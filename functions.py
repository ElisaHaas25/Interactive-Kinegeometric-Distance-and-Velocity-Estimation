import packages
from packages import *

# Source from R-file: call_velocity_prior.R

r_source = robjects.r['source']
r_source('./Rfiles/call_velocity_prior.R')

# Constants
kfac = 4.740471 # vel[km/s] = kfac * r[kpc] * propm[mas/yr]

def velpriorparams(dist,sfit):
    # Return the mean (meanTau) [1/(km/s)] and covariance (CovTau) [1/(km/s)^2]
    # of the velocity prior 
    
    # Returns the mean (meanTau) [1/(km/s)] and covariance (CovTau) [1/(km/s)^2]
    # of the velocity prior: 
    vramean_tau, vrasd_tau, vdecmean_tau, vdecsd_tau, cor_tau = robjects.r['eval.sfit'](sfit=sfit, r=dist*1e3) # output is in km/s, takes input in pc
    
    

    meanTau = np.array([vramean_tau,vdecmean_tau])
    CovTau = np.array([[vrasd_tau**2 ,vrasd_tau * vdecsd_tau *cor_tau], 
                       [vrasd_tau*vdecsd_tau*cor_tau,vdecsd_tau**2]])
    
    return meanTau , CovTau

def logdistpriordensity(dist, alpha,beta,rlen):
    # Return log (base 10) of the unnormalized distance prior [1/kpc]
    
    # convert rlen from pc to kpc 
    #rlen = rlen*1e-3

    #dist = np.where(dist > 0,dist,0)
    #prior = 1/(gamma((beta+1)/alpha)) * alpha/(rlen**(beta+1)) * dist**beta * np.exp(-((dist/rlen)**alpha))
    #logPrior = np.log10(prior)
    logPrior = -loggamma((beta+1)/alpha) + np.log(alpha) - (beta+1)*np.log(rlen) + beta*np.log(dist) - (dist/rlen)**alpha
    
    return 0.4342945*logPrior

def logparallaxlikelihood(dist, parallax, parallaxVar):
    # Return log (base 10) of the (normalized) likelhood of the parallax at this distance.
    # parallaxVar is the relevant element of the Cov3 matrix.
    # This is a simple 1D Gaussian density.
    #return 0.4342945* mvn.logpdf(parallax,mean=1/dist,cov=parallaxVar)
    return( 0.4342945*norm.logpdf(parallax,loc=1/dist,scale=np.sqrt(parallaxVar)) )
    
def loggeopostdensity(dist, parallax, parallaxVar, alpha,beta,rlen):
    # Return the log (base 10) of the unnormalized density of the geometric distance 
    # posterior [1/kpc]  
    # logparallaxlikelihood + logdistpriordensity
    
    #result = logparallaxlikelihood(dist=dist, parallax=parallax, parallaxVar=parallaxVar)\
    #+ logdistpriordensity(dist=dist,alpha=alpha,beta=beta,rlen=rlen)

    result = np.where(dist > 0, logparallaxlikelihood(dist=dist, parallax=parallax, 
                                                      parallaxVar=parallaxVar) + 
                       logdistpriordensity(dist=dist,alpha=alpha,beta=beta,rlen=rlen), -np.inf)

    return result

    
def logQfunc(dist, parallax, propm, Cov3, Cov2, kfac, meanTau, CovTau):
    """
    Implements Equation 8d
    In principle we only need the inverse of a covariance matrix to compute
    an unnormalized Gaussian, in which case we could just pass invCov2 and would
    not need Cov2. But scipy.stats.multivariate_normal doesn't seem to give this option.
    We could pass invCov2 to save inverting it for use in m_v and Cov_v, but it's only
    so just write it out explicitly.

    Parameters
    ----------
    dist:      distance [kpc] scalar
    parallax:  parallax [mas] scalar
    promp:     proper motion [mas/yr] 2-element vector
    Cov3:      full data astrometric covariance matrix [various] 3x3 matrix
    Cov2:      partial data covariance matrix [(mas/yr)^2] 2x2 matrix
    CovTau:    velocity prior covariance matrix [(km/s)^2] 2x2 matrix
    kfac:      constant

    Returns
    -------
    Log (base 10) density of Q(r,parallax,promp)
    """
    
    # Compute the inverse of Cov2 (equation 4b)
    # Compute X_mu (equation 4a)
    # Call velpriorparams() to get meanTau and CovTau
    # Invert CovTau
    # Compute logdensity of Q using scipy.stats.multivariate_normal.logpdf
    # See https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html
    # - don't forget to multiply by kfac*dist to get Q.
     
    Sigma_mu_w = np.array([Cov3[1,0],Cov3[2,0]])
    Sigma_w_w = Cov3[0,0]
    
    X_mu = Sigma_mu_w * Sigma_w_w**(-1) * (parallax-1/dist)
    
    
    
    logN = 0.4342945* mvn.logpdf(x=meanTau,mean=kfac*dist*(propm-X_mu),cov=kfac**2 * dist**2 * Cov2 + CovTau)
    
    return np.log10(kfac*dist) + logN
    
def logkinegeopostdensity(dist, parallax, parallaxVar ,propm, Cov3, Cov2, sfit, alpha,beta,rlen, kfac):
    # Return the log (base 10) of the (unnormalized) density of the kinegeometric 
    # distance posterior [1/kpc]  
    # logQfunc + loggeopostdensit
    
    meanTau, CovTau = velpriorparams(dist=dist, sfit=sfit) #
    
    #result = logQfunc(dist=dist, parallax=parallax, propm=propm, Cov3=Cov3, Cov2=Cov2, kfac=kfac, meanTau=meanTau, CovTau=CovTau)\
    #    + loggeopostdensity(dist=dist, parallax=parallax, parallaxVar=parallaxVar,alpha=alpha,beta=beta,rlen=rlen)
    
    result = np.where(dist > 0, logQfunc(dist=dist, parallax=parallax, propm=propm, Cov3=Cov3, Cov2=Cov2, kfac=kfac, meanTau=meanTau, CovTau=CovTau)
                      + loggeopostdensity(dist=dist, parallax=parallax, parallaxVar=parallaxVar,alpha=alpha,beta=beta,rlen=rlen),
                      -np.inf)
    return result

# Other functions only required for velocity posterior

def logvelpriordensity(vel, dist, sfit):
    # Return log (base 10) density of the velocity prior [1/(km/s)^2] using 
    # scipy.stats.multivariate_normal.logpdf
    # Call velpriorparams(dist, healpix)
    
    meanTau, CovTau = velpriorparams(dist=dist, sfit=sfit)
    N = mvn.pdf(x=vel,mean=meanTau,cov=CovTau)
    
    return np.log10(N)

def velpostparams(dist, parallax, Cov2, Cov3, propm, sfit, kfac): 
    # Return the mean (meanVel) [km/s] and covariance (CovVel) [km/s^2]
    # of the velocity posterior 
    
    meanTau, CovTau = velpriorparams(dist=dist, sfit=sfit)
    invCovTau = np.linalg.inv(CovTau)
    
    Sigma_mu_w = np.array([Cov3[1,0],Cov3[2,0]])
    Sigma_w_w = Cov3[0,0]
    
    X_mu = Sigma_mu_w * Sigma_w_w**(-1) * (parallax - (1/dist))
    mean2 = propm + X_mu
    
    A = np.linalg.inv(kfac**2 * dist**2 * Cov2)
    
    CovVel = np.linalg.inv(A + invCovTau)

    meanVel = np.matmul(CovVel , np.matmul(A*kfac*dist, propm - X_mu )+np.matmul(invCovTau,meanTau))
    
    
    return meanVel, CovVel

def plot_results(rSamp_kinegeo,rSamp_geo,rplotlo,rplothi,parallax, parallaxVar,
                 propm, Cov3, Cov2, sfit,alpha,beta,rlen, kfac, probs):
    
     # plot the results of the geometric and kinegeometric distance estimation                                                

    fig,ax = plt.subplots(3,1,figsize=(7,9), gridspec_kw={'height_ratios': [1,1,2]})
    
    fig.suptitle(f'Distance estimation ')
    
    r_plot = np.linspace(rplotlo,rplothi,200)
    
    # plot histogram of samples
    ax[0].hist(rSamp_kinegeo,bins=50,density=True,label='kinegeometric samples', color='black')
    ax[0].set_xlim(rplotlo,rplothi)
    ax[0].set_xlabel('distance [kpc]')
    ax[0].set_ylabel('density')
    
    # plot chains
    ax[1].plot(np.arange(0,len(rSamp_kinegeo),1),rSamp_kinegeo, color="black")
    ax[1].set_xlabel("iteration")
    ax[1].set_ylabel("kinegeo distance [kpc]")
    
    # plot kinegeo + geo posterior              
    logKinegeo_plot=[]
    logGeo_plot = []
    
    for i in r_plot: 
        logKinegeo_plot.append(logkinegeopostdensity(dist=i, parallax=parallax, parallaxVar=parallaxVar ,\
                                                     propm=propm, Cov3=Cov3, Cov2=Cov2, sfit=sfit,alpha=alpha,beta=beta,rlen=rlen, kfac=kfac))
        logGeo_plot.append(loggeopostdensity(dist=i, parallax=parallax, parallaxVar=parallaxVar, alpha=alpha,beta=beta,rlen=rlen))
        
    kinegeo_plot = 10**(np.array(logKinegeo_plot))         
    kinegeo_norm = integrate.trapezoid(kinegeo_plot,r_plot)
    
    geo_plot = 10**(np.array(logGeo_plot))
    geo_norm = integrate.trapezoid(geo_plot,r_plot)
    
    ax[2].plot(r_plot,kinegeo_plot/kinegeo_norm, label='kinegeometric posterior',color='orange')
    ax[0].plot(r_plot,kinegeo_plot/kinegeo_norm, label='kinegeometric posterior',color='orange')
    
    ax[2].plot(r_plot,geo_plot/geo_norm,label='geometric posterior',color='blue')
    
    ax[2].set_xlim(rplotlo,rplothi)
    ax[2].set_xlabel('distance [kpc]')
    ax[2].grid()
    
    # compute quantiles
    rQuant_kinegeo = np.quantile(rSamp_kinegeo,probs)
    rest_kinegeo = rQuant_kinegeo[0]
    rlo_kinegeo = rQuant_kinegeo[1]
    rhi_kinegeo = rQuant_kinegeo[2]
    rQuant_geo = np.quantile(rSamp_geo,probs)
    rest_geo = rQuant_geo[0]
    rlo_geo = rQuant_geo[1]
    rhi_geo = rQuant_geo[2]
    
    # overplot quantiles
    ax[2].axvline(rest_kinegeo,label ='$rmed_{kinegeo}$ (quantile 0.5)',color='orange')
    ax[2].axvline(rlo_kinegeo,linestyle='--',label ='$rlo_{kinegeo}$ (quantile 0.159)',color='orange')
    ax[2].axvline(rhi_kinegeo,linestyle='--',label ='$rhi_{kinegeo}$ (quantile 0.841)',color='orange')
    ax[2].axvline(rest_geo,label ='$rmed_{geo}$ (quantile 0.5)',color='blue')
    ax[2].axvline(rlo_geo,linestyle='--',label ='$rlo_{geo}$ (quantile 0.159)',color='blue')
    ax[2].axvline(rhi_geo,linestyle='--',label ='$rhi_{geo}$ (quantile 0.841)',color='blue')
    ax[2].legend(fontsize=7)
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    
    return fig


def plot_results_velocity(totVelsamp, probs):
    
    # plot results of velocity estimation: 
    
    fig, ax = plt.subplots(2,2)
    fig.suptitle(f'Velocity estimation')
    
    # Get results: mean and quantiles of v_ra and v_dec + correlation between array of v_ra and v_dec samples
    raVel = np.quantile(totVelsamp[:,0], probs)
    decVel = np.quantile(totVelsamp[:,1], probs)
    
    ax[0,0].hist(totVelsamp[:,0],histtype='step',color='k',bins=25)
    ax[0,0].sharex(ax[0,1])
    ax[0,0].axvline(raVel[0],0,1,color='k')
    ax[0,0].axvline(raVel[1],0,1,color='k',linestyle='--')
    ax[0,0].axvline(raVel[2],0,1,color='k',linestyle='--')
    
    ax[0,1].axis('off')
    ax[1,0].hist2d(totVelsamp[:,0], totVelsamp[:,1],cmap='Greys',bins=25)
    #ax[1,0].scatter(totVelsamp[:,0],totVelsamp[:,1],color='k',s=2)
    ax[1,0].set_xlabel('$v_{ra}$')
    ax[1,0].set_ylabel('$v_{dec}$')
    
    
    ax[1,1].hist(totVelsamp[:,1],orientation='horizontal', histtype='step',color='k',bins=25)
    ax[1,1].sharey(ax[0,1])
    ax[1,1].axhline(decVel[0],color='k')
    ax[1,1].axhline(decVel[1],color='k',linestyle='--')
    ax[1,1].axhline(decVel[2],color='k',linestyle='--')
    
    plt.setp(ax[0,0].get_xticklabels(), visible=False)
    plt.setp(ax[1,1].get_yticklabels(), visible=False);
    plt.tight_layout()
    
    return fig


# function to print the summary statistics for kinegeometric and geometric samples    
    
def print_summary_statistics(rInit,rStep,Nburnin,rSamp_kinegeo,rSamp_geo, totVelsamp, totMeanVel, n,thinfac,probs): 
    
    # print distance estimation statistics: 
    
    # compute quantiles of distance
    rQuant_kinegeo = np.quantile(rSamp_kinegeo,probs)
    rest_kinegeo = rQuant_kinegeo[0]
    rlo_kinegeo = rQuant_kinegeo[1]
    rhi_kinegeo = rQuant_kinegeo[2]
    rQuant_geo = np.quantile(rSamp_geo,probs)
    rest_geo = rQuant_geo[0]
    rlo_geo = rQuant_geo[1]
    rhi_geo = rQuant_geo[2]
    
    Nsamp=len(rSamp_kinegeo)
    
    print('\033[1m' + 'Distance estimation:' + '\033[0m')
    print('')
    print('MCMC initialization [pc]:', rInit)
    print('MCMC stepsize [pc]:',rStep)
    print('MCMC number of burn-in samples:',Nburnin)
    print('Thinning factor:', thinfac)
    print('MCMC number of retained iterations:',Nsamp)
    print('')
    print('Kinegeometric distance:')
    print('')
    print('estimated distance [pc] (quantile 0.5):',rest_kinegeo)
    print('lower distance limit [pc] (quantile 0.159):', rlo_kinegeo)
    print('upper distance limit [pc] (quantile 0.841):', rhi_kinegeo)
    print('')
    print('Geometric distance:')
    print('')
    print('estimated distance [pc] (quantile 0.5):',rest_geo)
    print('lower distance limit [pc] (quantile 0.159):', rlo_geo)
    print('upper distance limit [pc] (quantile 0.841):', rhi_geo)
    print('')
    
    #print velocity estimation statistics
    
    # Get results: mean and quantiles of v_ra and v_dec + correlation between array of v_ra and v_dec samples
    raVel = np.quantile(totVelsamp[:,0], probs)
    decVel = np.quantile(totVelsamp[:,1], probs)
    corrVel = np.corrcoef(totVelsamp[:,0], totVelsamp[:,1])
    
    # Expectation value of velocity
    E_v = 1/Nsamp * sum(totMeanVel)
    
    # Covariance between each velocity and distance
    r_mean = 1/Nsamp * sum(rSamp_kinegeo)
    Cov_rv_ra = 1/Nsamp * sum((rSamp_kinegeo-r_mean)*(totMeanVel[:,0]-E_v[0]))
    Cov_rv_dec = 1/Nsamp * sum((rSamp_kinegeo-r_mean)*(totMeanVel[:,1]-E_v[1]))
    
    print('\033[1m' + 'Velocity estimation:' + '\033[0m')
    print('')
    print('Number of velocity samples drawn for each of the MCMC distance samples: ',n)
    print('')
    print('Estimated velocities:')
    print('')
    print('velocity in ra:')
    print('estimated velocity [km/s] (quantile 0.5): ', raVel[0])
    print('lower velocity limit [km/s] (quantile 0.159): ',raVel[1])
    print('upper velocity limit [km/s](quantile 0.841): ',raVel[2])
    print('')
    print('velocity in dec:')
    print('estimated velocity [km/s] (quantile 0.5): ', decVel[0])
    print('lower velocity limit [km/s] (quantile 0.159): ',decVel[1])
    print('upper velocity limit [km/s](quantile 0.841): ',decVel[2])
    print('')
    print('Correlation between v_ra and v_dec: ', corrVel[1,0])
    print('')
    print('Covariance between each velocity and distance:' )
    print('Cov_rv_ra [kpc km/s]: ',Cov_rv_ra )
    print('Cov_rv_dec [kpc km/s]: ',Cov_rv_dec )

# function to print the data for a single source_id as input    
    
def print_data_sourceid(source_id,healpix,w,sd_w,wzp,parallax,mu_ra,mu_dec,sd_mu_ra,sd_mu_dec,corr_w_mu_ra,corr_w_mu_dec,corr_mu_ra_dec,alpha,beta,rlen):
    # print all the input data
    
    print('\033[1m' + f'Data for Gaia DR3 {source_id}: ' + '\033[0m')
    print('') 
    print('HEALpixel level 5:',healpix)
    print('alpha:',alpha)
    print('beta:', beta)
    print('rlen [kpc]:',rlen)
    print('')
    print('parallax [mas]',w)
    print('parallax error [mas]',sd_w)
    print('Zeropoint [mas]', wzp)
    print('Zeropoint-corrected parallax w[mas]',parallax)
    print('proper motion ra [mas/yr]:',mu_ra)
    print('proper motion dec [mas/yr]:',mu_dec)
    print('proper motion ra error [mas/yr]:',sd_mu_ra)
    print('proper motion dec error [mas/yr]:',sd_mu_dec)
    print('parallax - proper motion ra correlation: ',corr_w_mu_ra)
    print('parallax - proper motion dec correlation: ',corr_w_mu_dec)
    print('proper motion ra-dec correlation: ',corr_mu_ra_dec)
    print('')



    
# function that resolves a simbad name. If the name does not exist or there is no corresponding source_id, an error-string is returned. Only the source_id returned as string without the 'Gaia DR3' in front of it. The function at the top helps to find the Gaia DR3 name in the list of names from Simbad.

def extract_strings_with_word(string_list, word):
    result = [s for s in string_list if word in s]
    return result

def resolve_simbad_to_gaia(simbad_name):
    
    result_table = Simbad.query_objectids(simbad_name)
    
    if result_table is not None:
        
        gaia_dr3_source_id = extract_strings_with_word(result_table['ID'], 'Gaia DR3')
        
        if len(gaia_dr3_source_id) > 0:
            
            return re.findall(r'\d+',gaia_dr3_source_id[0])[1]
        
        else:
            
            return f'Error: Gaia DR3 source_id of {simbad_name} not found!'

    else:
        return f'Error: {simbad_name} not found in Simbad!'
from packages import *
from functions import *

def single_sourceid(sampler ,seed, source_id, 
                    parallax, parallax_error, alpha, beta, rlen, 
                    mu_ra, mu_dec, sd_mu_ra, sd_mu_dec, corr_w_mu_ra, corr_w_mu_dec, corr_mu_ra_dec, healpix,
                    Nwalker, a, Nsamp, thinfac, Nburnin, n, rows_prior_summary,probs): 
    
    np.random.seed(seed)
    
        
    if not source_id.isdigit():
    
        source_id = f'{resolve_simbad_to_gaia(source_id)}'
    
    if not source_id.isdigit(): 
    
        print(f"\033[1;31m {source_id} \033[0m") 
        
    else: 
    
        job = Gaia.launch_job("select "
                                            "source_id, parallax, parallax_error,\
                                            phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour,\
                                            ecl_lat, astrometric_params_solved,bp_rp,pmra,pmra_error,\
                                            pmdec, pmdec_error,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr "
                                            "from gaiadr3.gaia_source "
                                            "where source_id={}".format(int(source_id)))
        r = job.get_results()
        
        if len(r) == 0:
        
            print("\033[1;31m source_id does not exist! \033[0m")
        
        else: 
            w = float(r['parallax'])
            sd_w = float(r['parallax_error'])   
            mu_ra = float(r['pmra'])
            mu_dec = float(r['pmdec'])
            sd_mu_ra = float(r['pmra_error'])
            sd_mu_dec = float(r['pmdec_error'])
            corr_w_mu_ra = float(r['parallax_pmra_corr'])
            corr_w_mu_dec = float(r['parallax_pmdec_corr'])
            corr_mu_ra_dec = float(r['pmra_pmdec_corr'])
            phot_g_mean_mag = float(r['phot_g_mean_mag'])
            nu_eff_used_in_astrometry = float(r['nu_eff_used_in_astrometry'])
            pseudocolour = float(r['pseudocolour'])
            ecl_lat = float(r['ecl_lat'])
            astrometric_params_solved = float(r['astrometric_params_solved'])
            
            source_id = int(r['SOURCE_ID'])
            healpix = math.floor(source_id / (2**(35)*4**(12-5)) )
            
            # compute parallax zeropoint correction
            if astrometric_params_solved == 31 or astrometric_params_solved == 95:    
                if  phot_g_mean_mag == np.nan:
                    wzp = -0.017                    
                else:
                     wzp = zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)             
            else: 
                wzp = 0
            
            parallax = w - wzp
            parallaxVar = sd_w**2
            propm = np.array([mu_ra,mu_dec])
            Cov3 = np.array([[sd_w**2, sd_w*sd_mu_ra*corr_w_mu_ra, sd_w*sd_mu_dec*corr_w_mu_dec],
                             [sd_w*sd_mu_ra*corr_w_mu_ra, sd_mu_ra**2, sd_mu_ra*sd_mu_dec*corr_mu_ra_dec],
                             [sd_w*sd_mu_dec*corr_w_mu_dec, sd_mu_ra*sd_mu_dec*corr_mu_ra_dec, sd_mu_dec**2]])
            Sigma_w_w = Cov3[0,0]
            Sigma_mu_w = np.array([Cov3[1,0],Cov3[2,0]])
            Sigma_w_mu = np.transpose(Sigma_mu_w) #np.array([Cov3[0,1],Cov3[0,2]])
            Sigma_mu_mu = np.array([[Cov3[1,1],Cov3[1,2]],
                                   [Cov3[2,1],Cov3[2,2]]])
            Cov2 = Sigma_mu_mu - np.outer(Sigma_mu_w,np.dot(Sigma_w_w**(-1),Sigma_w_mu))
        
            alpha = float(rows_prior_summary[healpix][6])
            beta = float(rows_prior_summary[healpix][7])
            rlen = float(rows_prior_summary[healpix][5])*1e-3
            
            rlen_EDSD = float(rows_prior_summary[healpix][10])  # in pc, since only reqired function for mode_post3, which takes pc
            rInit = float(mode_post3(w=1e-3*parallax,wsd=1e-3*sd_w,rlen = rlen_EDSD,retall = False))*1e-3 # in kpc
            rStep = 0.75*rInit*min(1/3, abs(sd_w/parallax))
            
            
            
            # load the velocity prior model (sfit) correspoding to the respective healpixel 
            
            # Activate R's base package
            base = importr('base')
            
            # Create a new R environment
            tempEnv = robjects.Environment()
            
            # Load the R object from the URL
            url = "https://www2.mpia-hd.mpg.de/homes/calj/gedr3_distances_velocities/velocity_prior_fits/models/"
            
            # Load the R object into tempEnv
            robjects.r['load'](robjects.r.url(url + str(healpix) + ".Robj"), envir=tempEnv)
            
            # Extract sfit from tempEnv
            sfit = tempEnv['sfit']
            
            # Distance estimation
            
            walker_init = np.random.uniform((1-a)*rInit, (1+a)*rInit, (Nwalker, 1))
            
            if sampler=='metropolis': 
                
                # geometric samples
                
                samp_geo = metrop(func=loggeopostdensity, thetaInit=rInit, Nburnin=Nburnin, Nsamp=Nsamp, sampleCov=rStep**2, 
                              parallax=parallax, parallaxVar=parallaxVar, alpha=alpha, beta=beta, rlen=rlen, seed=seed)
                
                rSamp_geo = samp_geo[:,1]
                rSamp_geo = rSamp_geo[::thinfac]
                
                # kinegeometric samples
                
                samp_kinegeo = metrop(func=logkinegeopostdensity, thetaInit=rInit, Nburnin=Nburnin, Nsamp=Nsamp, sampleCov=rStep**2, 
                                  parallax=parallax, parallaxVar=parallaxVar, propm=propm, Cov3=Cov3, Cov2=Cov2, 
                                  sfit=sfit, alpha=alpha, beta=beta, rlen=rlen, kfac=kfac, seed=seed)
                
                rSamp_kinegeo = samp_kinegeo[:,1]
                rSamp_kinegeo = rSamp_kinegeo[::thinfac]
                
            if sampler == 'emcee': 
                
                # geometric samples
                
                sampler_geo = emcee.EnsembleSampler(nwalkers=Nwalker, 
                                            ndim=1, 
                                            log_prob_fn=loggeopostdensity, 
                                            args=[parallax, parallaxVar, alpha,beta,rlen])
                sampler_geo.run_mcmc(walker_init, Nsamp+Nburnin)
                samp_geo = sampler_geo.get_chain(flat=True, discard = Nburnin) 
                rSamp_geo = samp_geo.flatten()
                rSamp_geo = rSamp_geo[::thinfac]
                
                # kinegeometric samples
                sampler_kinegeo = emcee.EnsembleSampler(nwalkers=Nwalker,
                                                        ndim=1,
                                                        log_prob_fn=logkinegeopostdensity,
                                                        args=[parallax, parallaxVar ,propm, Cov3, Cov2, sfit,alpha,beta,rlen, kfac])
                
                sampler_kinegeo.run_mcmc(walker_init,Nsamp+Nburnin)
                samp_kinegeo = sampler_kinegeo.get_chain(flat=True,discard = Nburnin) 
                rSamp_kinegeo = samp_kinegeo.flatten()
                rSamp_kinegeo = rSamp_kinegeo[::thinfac]
           
            
            # Velocity sampling: 
            
            # sample n velocities and take the mean for each kinegeo distance sample

            totVelsamp = [] # mean of velocity samples for all distances
            totMeanVel = [] # mean of velocity posterior for all distances 
            
            for i in range(len(rSamp_kinegeo)): 
                meanVel, CovVel = velpostparams(dist=rSamp_kinegeo[i], parallax=parallax, Cov2=Cov2, Cov3=Cov3, propm=propm, sfit=sfit,
                                                kfac=kfac)
                velsamp = np.random.multivariate_normal(mean = meanVel, cov = CovVel,size = n)
                velsamp_mean = np.mean(velsamp, axis=0)
                totVelsamp.append(velsamp_mean)
                totMeanVel.append(meanVel)
            totVelsamp = np.array(totVelsamp)
            totMeanVel = np.array(totMeanVel)
            
            #print data 
            print_data_sourceid(source_id=source_id,healpix=healpix,w=w,sd_w=sd_w,
                                wzp=wzp,parallax=parallax,mu_ra=mu_ra,mu_dec=mu_dec,sd_mu_ra=sd_mu_ra,
                                sd_mu_dec=sd_mu_dec,corr_w_mu_ra=corr_w_mu_ra,corr_w_mu_dec=corr_w_mu_dec,corr_mu_ra_dec=corr_mu_ra_dec,
                               alpha=alpha, beta=beta, rlen=rlen)
            
            #plot results
            rplotlo = min(np.percentile(rSamp_kinegeo,5),np.percentile(rSamp_geo,5))
            rplothi = max(np.percentile(rSamp_kinegeo,95),np.percentile(rSamp_geo,95))
            
            plot_results(rSamp_kinegeo=rSamp_kinegeo,rSamp_geo=rSamp_geo,
                         rplotlo=rplotlo,rplothi=rplothi,parallax=parallax, parallaxVar=parallaxVar, 
                         propm=propm, Cov3=Cov3, Cov2=Cov2, sfit=sfit,alpha=alpha,beta=beta,rlen=rlen, kfac=kfac,probs=probs) 
            
            plot_results_velocity(totVelsamp=totVelsamp,probs=probs)
            
            #print summary statistics
            print_summary_statistics(sampler=sampler, walker_init=walker_init,
                                     rInit=rInit,rStep=rStep,Nburnin=Nburnin,rSamp_kinegeo=rSamp_kinegeo,rSamp_geo=rSamp_geo,
                                     totVelsamp=totVelsamp, totMeanVel=totMeanVel,n=n,thinfac=thinfac,probs=probs)
            plt.show()
                
                
def single_owndata(sampler,seed, source_id, 
                   parallax, parallax_error, alpha, beta, rlen, 
                   mu_ra, mu_dec, sd_mu_ra, sd_mu_dec, corr_w_mu_ra, corr_w_mu_dec, corr_mu_ra_dec, healpix,         
                   Nwalker,a ,rInit ,rStep, Nsamp, thinfac, Nburnin, n,rows_prior_summary, probs):       
    
    np.random.seed(seed)
    
  
        
    propm = np.array([mu_ra,mu_dec])
    sd_w = parallax_error
    parallaxVar = sd_w **2
    Cov3 = np.array([[sd_w**2, sd_w*sd_mu_ra*corr_w_mu_ra, sd_w*sd_mu_dec*corr_w_mu_dec],
                 [sd_w*sd_mu_ra*corr_w_mu_ra, sd_mu_ra**2, sd_mu_ra*sd_mu_dec*corr_mu_ra_dec],
                 [sd_w*sd_mu_dec*corr_w_mu_dec, sd_mu_ra*sd_mu_dec*corr_mu_ra_dec, sd_mu_dec**2]])
    Sigma_w_w = Cov3[0,0]
    Sigma_mu_w = np.array([Cov3[1,0],Cov3[2,0]])
    Sigma_w_mu = np.transpose(Sigma_mu_w) #np.array([Cov3[0,1],Cov3[0,2]])
    Sigma_mu_mu = np.array([[Cov3[1,1],Cov3[1,2]],
                           [Cov3[2,1],Cov3[2,2]]])
    Cov2 = Sigma_mu_mu - np.outer(Sigma_mu_w,np.dot(Sigma_w_w**(-1),Sigma_w_mu))
        
    # load the velocity prior model (sfit) correspoding to the respective healpixel 
    
    # Activate R's base package
    base = importr('base')
    
    # Create a new R environment
    tempEnv = robjects.Environment()
    
    # Load the R object from the URL
    url = "https://www2.mpia-hd.mpg.de/homes/calj/gedr3_distances_velocities/velocity_prior_fits/models/"
    
    # Load the R object into tempEnv
    robjects.r['load'](robjects.r.url(url + str(healpix) + ".Robj"), envir=tempEnv)
    
    # Extract sfit from tempEnv
    sfit = tempEnv['sfit']
    
    walker_init = np.random.uniform((1-a)*rInit, (1+a)*rInit, (Nwalker, 1))
    
    if sampler=='metropolis': 
        
        # geometric samples
        
        samp_geo = metrop(func=loggeopostdensity, thetaInit=rInit, Nburnin=Nburnin, Nsamp=Nsamp, sampleCov=rStep**2, 
                      parallax=parallax, parallaxVar=parallaxVar, alpha=alpha, beta=beta, rlen=rlen, seed=seed)
        
        rSamp_geo = samp_geo[:,1]
        rSamp_geo = rSamp_geo[::thinfac]
        
        # kinegeometric samples
        
        samp_kinegeo = metrop(func=logkinegeopostdensity, thetaInit=rInit, Nburnin=Nburnin, Nsamp=Nsamp, sampleCov=rStep**2, 
                          parallax=parallax, parallaxVar=parallaxVar, propm=propm, Cov3=Cov3, Cov2=Cov2, 
                          sfit=sfit, alpha=alpha, beta=beta, rlen=rlen, kfac=kfac, seed=seed)
        
        rSamp_kinegeo = samp_kinegeo[:,1]
        rSamp_kinegeo = rSamp_kinegeo[::thinfac]
        
    if sampler == 'emcee': 
        
        # geometric samples
        
        sampler_geo = emcee.EnsembleSampler(nwalkers=Nwalker, 
                                    ndim=1, 
                                    log_prob_fn=loggeopostdensity, 
                                    args=[parallax, parallaxVar, alpha,beta,rlen])
        sampler_geo.run_mcmc(walker_init, Nsamp+Nburnin)
        samp_geo = sampler_geo.get_chain(flat=True, discard = Nburnin) 
        rSamp_geo = samp_geo.flatten()
        rSamp_geo = rSamp_geo[::thinfac]
        
        # kinegeometric samples
        
        sampler_kinegeo = emcee.EnsembleSampler(nwalkers=Nwalker,
                                                ndim=1,
                                                log_prob_fn=logkinegeopostdensity,
                                                args=[parallax, parallaxVar ,propm, Cov3, Cov2, sfit,alpha,beta,rlen, kfac])
        
        sampler_kinegeo.run_mcmc(walker_init,Nsamp+Nburnin)
        samp_kinegeo = sampler_kinegeo.get_chain(flat=True,discard = Nburnin) 
        rSamp_kinegeo = samp_kinegeo.flatten()
        rSamp_kinegeo = rSamp_kinegeo[::thinfac]
    
    # Velocity sampling: 
            
            # sample n velocities and take the mean for each kinegeo distance sample

    totVelsamp = [] # mean of velocity samples for all distances
    totMeanVel = [] # mean of velocity posterior for all distances 
    
    for i in range(len(rSamp_kinegeo)): 
        meanVel, CovVel = velpostparams(dist=rSamp_kinegeo[i], parallax=parallax, Cov2=Cov2, Cov3=Cov3, propm=propm, sfit=sfit, kfac=kfac)
        velsamp = np.random.multivariate_normal(mean = meanVel, cov = CovVel,size = n)
        velsamp_mean = np.mean(velsamp, axis=0)
        totVelsamp.append(velsamp_mean)
        totMeanVel.append(meanVel)
    totVelsamp = np.array(totVelsamp)
    totMeanVel = np.array(totMeanVel)
      
   #display results 
                
    rplotlo = min(np.percentile(rSamp_kinegeo,5),np.percentile(rSamp_geo,5))
    rplothi = max(np.percentile(rSamp_kinegeo,95),np.percentile(rSamp_geo,95))
    
    plot_results(rSamp_kinegeo=rSamp_kinegeo,rSamp_geo=rSamp_geo,
                 rplotlo=rplotlo,rplothi=rplothi,parallax=parallax, parallaxVar=parallaxVar, 
                 propm=propm, Cov3=Cov3, Cov2=Cov2, sfit=sfit,alpha=alpha,beta=beta,rlen=rlen, kfac=kfac,probs=probs) 
    plot_results_velocity(totVelsamp=totVelsamp,probs=probs)
    
    print_summary_statistics(sampler=sampler,walker_init=walker_init,rInit=rInit,rStep=rStep,Nburnin=Nburnin,rSamp_kinegeo=rSamp_kinegeo,rSamp_geo=rSamp_geo, 
                             totVelsamp=totVelsamp, totMeanVel=totMeanVel, n=n,thinfac=thinfac,probs=probs)
    
    plt.show()
    

def multiple_sorceid(filename_in,filename_out,create_pdf, sampler, Nsamp,Nwalker, a ,Nburnin,thinfac,n,seed,rows_prior_summary,Nmax,probs): 
    
    #read in comparison data to obtain array containing source_ids 
    
    source_ids = np.genfromtxt(f'./data/{filename_in}.csv', delimiter=',', skip_header=1,  dtype='int64', usecols=0) 
    print('Number of source_ids read from file: ',len(source_ids))
    
    #extract data from Gaia 
    source_ids_string = ", ".join(map(str, source_ids))
    job = Gaia.launch_job("select "
                                        "source_id, parallax, parallax_error,phot_g_mean_mag,\
                                        nu_eff_used_in_astrometry, pseudocolour,\
                                        ecl_lat, astrometric_params_solved,bp_rp,\
                                        pmra,pmra_error,pmdec, pmdec_error,parallax_pmra_corr,\
                                        parallax_pmdec_corr,pmra_pmdec_corr "
                                        "from gaiadr3.gaia_source "
                                        "where source_id in ({})".format(source_ids_string))
    r = job.get_results()
    
    print('Number of sources found on GACS: ', len(r))
    
    if len(r) > Nmax: 
        print(f"\033[93m Warning: processing more than {Nmax} sources.\033[0m")
        
    # samples for each source: 
    
    rSamples_kinegeo = []
    rSamples_geo = []
    
    
    
    # distance statistics for each source
    
    rMedKinegeo_all = []
    rLoKinegeo_all = []
    rHiKinegeo_all = []
    
    rMedGeo_all = []
    rLoGeo_all = []
    rHiGeo_all = []
    
    
    # velocity statistics for each source:
    
    raVel_all = []
    decVel_all = []
    corrVel_all = []
    Cov_rv_ra_all = []
    Cov_rv_dec_all = []

    totVelsamp_all = []

    # plots for each source: 
    
    plots_dist = []
    plots_vel = []
        
    for i in range(len(source_ids[:])): 
    
        w = float(r['parallax'][i])
        sd_w = float(r['parallax_error'][i])   
        
        mu_ra = float(r['pmra'][i])
        mu_dec = float(r['pmdec'][i])
        sd_mu_ra = float(r['pmra_error'][i])
        sd_mu_dec = float(r['pmdec_error'][i])
        corr_w_mu_ra = float(r['parallax_pmra_corr'][i])
        corr_w_mu_dec = float(r['parallax_pmdec_corr'][i])
        corr_mu_ra_dec = float(r['pmra_pmdec_corr'][i])
        
        phot_g_mean_mag = float(r['phot_g_mean_mag'][i])
        nu_eff_used_in_astrometry = float(r['nu_eff_used_in_astrometry'][i])
        pseudocolour = float(r['pseudocolour'][i])
        ecl_lat = float(r['ecl_lat'][i])
        astrometric_params_solved = float(r['astrometric_params_solved'][i])
        
        
        source_id = int(r['SOURCE_ID'][i])
        print(f'Processing star {i+1} of {len(r)}: Gaia DR3 {source_id}')
        
        healpix = math.floor(source_id / (2**(35)*4**(12-5)) )
        
        # zeropoint correction
        if astrometric_params_solved == 31 or astrometric_params_solved == 95:    
            if  phot_g_mean_mag == np.nan:
                wzp = -0.017                    
            else:
                 wzp = zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)             
        else: 
            wzp = 0
        
        parallax = w - wzp
        parallaxVar = sd_w**2
        propm = np.array([mu_ra,mu_dec])
        
        Cov3 = np.array([[sd_w**2, sd_w*sd_mu_ra*corr_w_mu_ra, sd_w*sd_mu_dec*corr_w_mu_dec],
                         [sd_w*sd_mu_ra*corr_w_mu_ra, sd_mu_ra**2, sd_mu_ra*sd_mu_dec*corr_mu_ra_dec],
                         [sd_w*sd_mu_dec*corr_w_mu_dec, sd_mu_ra*sd_mu_dec*corr_mu_ra_dec, sd_mu_dec**2]])
        
        
        Sigma_w_w = Cov3[0,0]
        Sigma_mu_w = np.array([Cov3[1,0],Cov3[2,0]])
        Sigma_w_mu = np.transpose(Sigma_mu_w) #np.array([Cov3[0,1],Cov3[0,2]])
        Sigma_mu_mu = np.array([[Cov3[1,1],Cov3[1,2]],
                               [Cov3[2,1],Cov3[2,2]]])
        
        
        
        Cov2 = Sigma_mu_mu - np.outer(Sigma_mu_w, Sigma_w_w**(-1)*Sigma_w_mu)
        
        
        alpha = float(rows_prior_summary[healpix][6])
        beta = float(rows_prior_summary[healpix][7])
        rlen = float(rows_prior_summary[healpix][5])*1e-3 # in kpc
        rlen_EDSD = float(rows_prior_summary[healpix][10])  # in pc, since only reqired function for mode_post3, which takes pc
        
        rInit = float(mode_post3(w=1e-3*parallax,wsd=1e-3*sd_w,rlen = rlen_EDSD,retall = False))*1e-3 # in kpc
        rStep = 0.75*rInit*min(1/3, abs(sd_w/parallax))
        
        # load the velocity prior model (sfit) correspoding to the respective healpixel 
    
        # Activate R's base package
        base = importr('base')
        
        # Create a new R environment
        tempEnv = robjects.Environment()
        
        # Load the R object from the URL
        url = "https://www2.mpia-hd.mpg.de/homes/calj/gedr3_distances_velocities/velocity_prior_fits/models/"
        
        # Load the R object into tempEnv
        robjects.r['load'](robjects.r.url(url + str(healpix) + ".Robj"), envir=tempEnv)
        
        # Extract sfit from tempEnv
        sfit = tempEnv['sfit']
        
        if sampler == 'metropolis':
        # Uses metropolis algorithm from metropolis.py
        
            #kinegeometric samples
            
            samp_kinegeo = metrop(func=logkinegeopostdensity ,thetaInit= rInit ,Nburnin=Nburnin ,Nsamp=Nsamp,
                                  sampleCov=rStep**2 ,parallax=parallax, seed=seed,parallaxVar = parallaxVar,propm=propm, Cov3=Cov3,
                                  Cov2=Cov2, sfit=sfit, alpha=alpha,beta=beta,rlen=rlen,kfac=kfac)
            rSamp_kinegeo = samp_kinegeo[:,1]
            rSamp_kinegeo = rSamp_kinegeo[::thinfac]
            
            rSamples_kinegeo.append(rSamp_kinegeo)
            
            #geometric samples
            samp_geo = metrop(func=loggeopostdensity ,thetaInit= rInit ,Nburnin=Nburnin ,Nsamp=Nsamp,sampleCov=rStep**2 ,seed=seed,
                              parallax=parallax, parallaxVar = parallaxVar,alpha=alpha,beta=beta,rlen=rlen)
            rSamp_geo = samp_geo[:,1]
            rSamp_geo = rSamp_geo[::thinfac]
            rSamples_geo.append(rSamp_geo)
        
        if sampler == 'emcee': 
            
            walker_init = np.random.uniform((1-a)*rInit, (1+a)*rInit, (Nwalker, 1))
            
            # geometric samples
            sampler_geo = emcee.EnsembleSampler(nwalkers=Nwalker, 
                                        ndim=1, 
                                        log_prob_fn=loggeopostdensity, 
                                        args=[parallax, parallaxVar, alpha,beta,rlen])
            sampler_geo.run_mcmc(walker_init, Nsamp+Nburnin)
            samp_geo = sampler_geo.get_chain(flat=True, discard = Nburnin) 
            rSamp_geo = samp_geo.flatten()
            rSamp_geo = rSamp_geo[::thinfac]
            
            # kinegeometric samples
            sampler_kinegeo = emcee.EnsembleSampler(nwalkers=Nwalker,
                                                    ndim=1,
                                                    log_prob_fn=logkinegeopostdensity,
                                                    args=[parallax, parallaxVar ,propm, Cov3, Cov2, sfit,alpha,beta,rlen, kfac])
            
            sampler_kinegeo.run_mcmc(walker_init,Nsamp+Nburnin)
            samp_kinegeo = sampler_kinegeo.get_chain(flat=True,discard = Nburnin) 
            rSamp_kinegeo = samp_kinegeo.flatten()
            rSamp_kinegeo = rSamp_kinegeo[::thinfac]
            
            
        #compute distance statistics:
                 
        rQuantKinegeo = np.quantile(rSamp_kinegeo,probs)
        
        rMedKinegeo = rQuantKinegeo[0]
        rLoKinegeo = rQuantKinegeo[1]
        rHiKinegeo= rQuantKinegeo[2]
        
        rMedKinegeo_all.append(rMedKinegeo)
        rLoKinegeo_all.append(rLoKinegeo)
        rHiKinegeo_all.append(rHiKinegeo)
        
        rQuantGeo = np.quantile(rSamp_geo,probs)
        
        rMedGeo = rQuantGeo[0]
        rLoGeo = rQuantGeo[1]
        rHiGeo= rQuantGeo[2]
        
        rMedGeo_all.append(rMedGeo)
        rLoGeo_all.append(rLoGeo)
        rHiGeo_all.append(rHiGeo)
        
        
        # compute velocity statisitics: 
            
        # loop over all samples of a source:
        
        totVelsamp = [] #mean of velocity samples for all distances
        totMeanVel = [] #mean of velocity posterior for all distances 
    
        for i in range(len(rSamp_kinegeo)): 
            
            meanVel, CovVel = velpostparams(dist=rSamp_kinegeo[i], parallax=parallax, Cov2=Cov2, Cov3=Cov3, propm=propm, sfit=sfit,
                                            kfac=kfac)
            velsamp = np.random.multivariate_normal(mean = meanVel, cov = CovVel,size = n)
            velsamp_mean = np.mean(velsamp, axis=0)
            totVelsamp.append(velsamp_mean)
            totMeanVel.append(meanVel)
            
        totVelsamp_all.append(totVelsamp)
        totVelsamp = np.array(totVelsamp)
        totMeanVel = np.array(totMeanVel)
        
        #get results: mean and quantiles of v_ra and v_dec + correlation between array of v_ra and v_dec samples
            
        raVel = np.quantile(totVelsamp[:,0], probs)
        decVel = np.quantile(totVelsamp[:,1], probs)
        corrVel = np.corrcoef(totVelsamp[:,0], totVelsamp[:,1])
        
        # Expectation value of velocity
            
        E_v = 1/Nsamp * sum(totMeanVel)
        
        #covariance between each velocity and distance
        
        r_mean = 1/Nsamp * sum(rSamp_kinegeo)
        
        Cov_rv_ra = 1/Nsamp * sum((rSamp_kinegeo-r_mean)*(totMeanVel[:,0]-E_v[0]))
        Cov_rv_dec = 1/Nsamp * sum((rSamp_kinegeo-r_mean)*(totMeanVel[:,1]-E_v[1]))
            
    
        raVel_all.append(raVel)
        decVel_all.append(decVel)
        corrVel_all.append(corrVel)
        Cov_rv_ra_all.append(Cov_rv_ra)
        Cov_rv_dec_all.append(Cov_rv_dec)
    
        robjects.r['closeAllConnections']()
            
        # save to csv file; all samples for one source_id are in one row, same order as input 
        
        np.savetxt(f'./results/{filename_out}_samples_kinegeo.csv', rSamples_kinegeo, delimiter=",")
        np.savetxt(f'./results/{filename_out}_samples_geo.csv', rSamples_geo, delimiter=",")
        
        # create plots to save in pdf file and append to the lists of plots for each source: 
        
        if create_pdf == True:
            
            rplotlo = min(np.percentile(rSamp_kinegeo,5),np.percentile(rSamp_geo,5))
            rplothi = max(np.percentile(rSamp_kinegeo,95),np.percentile(rSamp_geo,95))

            #rplotlo = 0.2*min(rSamp_kinegeo)
            #rplothi = 1.2*max(rSamp_kinegeo)
            
            fig_dist = plot_results(rSamp_kinegeo=rSamp_kinegeo,rSamp_geo=rSamp_geo,
                                    rplotlo=rplotlo,rplothi=rplothi,parallax=parallax, parallaxVar=parallaxVar, 
                                    propm=propm, Cov3=Cov3, Cov2=Cov2, sfit=sfit,alpha=alpha,beta=beta,rlen=rlen, kfac=kfac,probs=probs)     
            fig_dist.suptitle(f'Distance estimation of Gaia DR3 {source_id}')
            
            fig_vel = plot_results_velocity(totVelsamp=totVelsamp,probs=probs)    
            
            fig_vel.suptitle(f'Velocity estimation of Gaia DR3 {source_id}')
            
            plots_dist.append(fig_dist)
            plots_vel.append(fig_vel)
        
    # save all summary statistics to a csv file

    header = ['source_id',\
              'rMedGeo', 'rLoGeo', 'rHiGeo',\
              'rMedKinogeo', 'rLoKinogeo', 'rHiKinogeo',\
             'vRaMedKinogeo', 'vRaLoKinogeo', 'vRaHiKinogeo',\
             'vDecMedKinogeo', 'vDecLoKinogeo', 'vDecHiKinogeo',\
             'rvraCorrKinogeo', 'rvdecCorrKinogeo', 'vravdecCorrKinogeo']
    
    with open(f'./results/{filename_out}_summary_statistics.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for i in range(len(source_ids[:])): 
            writer.writerow([source_ids[i],\
                             rMedGeo_all[i], rLoGeo_all[i], rHiGeo_all[i],\
                             rMedKinegeo_all[i], rLoKinegeo_all[i], rHiKinegeo_all[i],\
                             raVel_all[i][0],raVel_all[i][1],raVel_all[i][2],\
                             decVel_all[i][0],decVel_all[i][1],decVel_all[i][2],\
                             Cov_rv_ra_all[i],Cov_rv_dec_all[i],corrVel_all[i][0,1]])  
            
    # save all plots to a .pdf file:
    
    if create_pdf == True:
    
        with PdfPages(f'./results/{filename_out}.pdf') as pdf:
            
            for i in range(len(plots_dist)):
                pdf.savefig(plots_dist[i])
                pdf.savefig(plots_vel[i])
                        
def interactive_velocity_function(data, sampler, seed, source_id,
                                  parallax, parallax_error, alpha, beta, rlen, 
                                  mu_ra, mu_dec, sd_mu_ra, sd_mu_dec, corr_w_mu_ra, corr_w_mu_dec, corr_mu_ra_dec, healpix,
                                  Nwalker,a, rInit, rStep, Nsamp, thinfac, Nburnin, n,
                                  filename_in, filename_out, create_pdf, rows_prior_summary, Nmax, probs):
    if Nsamp < 0: 
        
        print("\033[91mError: number of posterior samples must not be negative.\033[0m")
        
    elif Nburnin < 0: 
        
        print("\033[91mError: number of burn-in samples must not be negative.\033[0m")
    
    elif thinfac <= 0: 
        
        print("\033[91mError: thinning factor must not be less or equal 0.\033[0m")
   
    elif n < 0: 
        
        print("\033[91mError: Number of velocity samples drawn for each of the MCMC distance samples must nut be negative.\033[0m")
    
    else: 
        message= 'Processing...'
        print(message, end='\r')
        
        if data == 'single source, own data': 
            
            single_owndata(sampler=sampler,seed=seed, source_id=source_id,
                           parallax=parallax, parallax_error=parallax_error, alpha=alpha, beta=beta, rlen=rlen,
                           mu_ra=mu_ra, mu_dec=mu_dec, sd_mu_ra=sd_mu_ra, sd_mu_dec=sd_mu_dec,                             
                           corr_w_mu_ra=corr_w_mu_ra, corr_w_mu_dec=corr_w_mu_dec, corr_mu_ra_dec=corr_mu_ra_dec,                         
                           healpix=healpix,
                           Nwalker=Nwalker,a=a, rInit=rInit, rStep=rStep, Nsamp=Nsamp, thinfac=thinfac, Nburnin=Nburnin,n=n,
                           rows_prior_summary=rows_prior_summary, probs=probs)
                               
            
        if data == 'single source, source_id': 
            
            single_sourceid(sampler=sampler,seed=seed, source_id=source_id,
                            parallax=parallax, parallax_error=parallax_error, alpha=alpha, beta=beta, rlen=rlen,
                            mu_ra=mu_ra, mu_dec=mu_dec, sd_mu_ra=sd_mu_ra, sd_mu_dec=sd_mu_dec,                             
                            corr_w_mu_ra=corr_w_mu_ra, corr_w_mu_dec=corr_w_mu_dec, corr_mu_ra_dec=corr_mu_ra_dec,                         
                            healpix=healpix,
                            Nwalker=Nwalker,a=a, Nsamp=Nsamp, thinfac=thinfac, Nburnin=Nburnin,n=n,
                            rows_prior_summary=rows_prior_summary, probs=probs)                           
                                                               
        if data =='multiple source_ids in .csv file':
            
            multiple_sorceid(filename_in=filename_in,filename_out=filename_out,create_pdf=create_pdf,
                             Nwalker=Nwalker,a=a, sampler=sampler, Nsamp=Nsamp, Nburnin=Nburnin,
                             thinfac=thinfac, n=n,seed=seed,
                             rows_prior_summary=rows_prior_summary,Nmax=Nmax, probs=probs)
            
        print(" " * len(message), end='\r')
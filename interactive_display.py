from packages import *
from functions import *
from interactive_functions import *
from configurations import *
rows_prior_summary = np.loadtxt('prior_summary.csv', delimiter=',',skiprows=1)

#%matplotlib inline

seed = widgets.IntText(
    value=123,
    description='seed',
    disabled=False,
    style = {'description_width': 'initial'}
)

#button to select the sampler

sampler = widgets.RadioButtons(
    options=['metropolis', 'emcee'],
    description='Sampler:',
    disabled=False
)

#field to insert source_id

source_id = widgets.Text(
    value='3490289711213205632',
    description='source_id/name',
    disabled=False,
    style = {'description_width': 'initial'}
)

description_text = widgets.Label(value='Gaia DR3 source_id (without ‘Gaia DR3’) or name to be resolved at Simbad.')

#field to insert name of file containing source_ids or (name, parallax, parallax_error, ra, dec)

filename_in = widgets.Text(
    value='velocity_test',
    description='Input .csv-file',
    disabled=False
)

#field to name the output file

filename_out = widgets.Text(
    value='velocity_results',
    description='Output .pdf/.csv-file',
    disabled=False,
    style = {'description_width': 'initial'}
)

# Buttons to select the mode

data = widgets.Select(
    options=['single source, own data', 'single source, source_id', 'multiple source_ids in .csv file'],
    value='single source, own data',
    # rows=10,
    description='Mode:',
    layout={'width': 'max-content'},
    disabled=False
)


# rInit slider

rInit0 = widgets.FloatSlider(
    value=1,
    min=0,
    max=2,
    step=0.001,
    description='starting value [kpc]:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.4f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax metrop-start_slider

rInit = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((rInit0, 'value'), (rInit, 'value'))

#stepsize slider

rStep0 = widgets.FloatSlider(
    value=0.25,
    min=0,
    max=0.5,
    step=0.001,
    description='stepsize (only used for Metropolis) [kpc]:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.3f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to rStep_slider

rStep = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((rStep0, 'value'), (rStep, 'value'))

#Nsamp slider

Nsamp0 = widgets.IntSlider(
    value = 5000,
    min=0,
    max=10000,
    step=1,
    description='number of posterior samples after removing burnins:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax Nsamp_slider

Nsamp = widgets.IntText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((Nsamp0, 'value'), (Nsamp, 'value'))

#n slider

n0 = widgets.IntSlider(
    value = 1,
    min=1,
    max=10,
    step=1,
    description='Number of velocity samples drawn for each of the MCMC distance samples:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

n = widgets.IntText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((n0, 'value'), (n, 'value'))


#tninning slider

thinfac0 = widgets.IntSlider(
    value = 10,
    min=0,
    max=1000,
    step=1,
    description='thinning factor:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax Nsamp_slider

thinfac = widgets.IntText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((thinfac0, 'value'), (thinfac, 'value'))

# walker_init

#Nburnin slider

Nburnin0 = widgets.IntSlider(
    value = 500,
    min=0,
    max=1000,
    step=1,
    description='number of burn-in samples:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to Nburnin slider

Nburnin = widgets.IntText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((Nburnin0, 'value'), (Nburnin, 'value'))

#Nwalker slider

Nwalker0 = widgets.IntSlider(
    value = 5,
    min=0,
    max=20,
    step=1,
    description='Number of walkers:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to Nwalker slider

Nwalker = widgets.IntText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((Nwalker0, 'value'), (Nwalker, 'value'))

#a slider

a0 = widgets.FloatSlider(
    value = 0.01,
    min=0,
    max=1,
    step=0.005,
    description='a:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax slider

a = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((a0, 'value'), (a, 'value'))


#parallax slider

parallax0 = widgets.FloatSlider(
    value = 1,
    min=-0.5,
    max=3,
    step=0.01,
    description='parallax [mas]:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax slider

parallax = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((parallax0, 'value'), (parallax, 'value'))

#parallax_error slider

parallax_error0 = widgets.FloatSlider(
    value = 0.5,
    min=0.001,
    max=3,
    step=0.001,
    description='parallax error [mas]:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax_error slider

parallax_error = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((parallax_error0, 'value'), (parallax_error, 'value'))

#rlen slider

rlen0 = widgets.FloatSlider(
    value= 5,
    min=0,
    max=10,
    step=0.001,
    description='Length scale [kpc]:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax rlen-slider

rlen = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((rlen0, 'value'), (rlen, 'value'))

# alpha-slider

alpha0 = widgets.FloatSlider(
    value=1,
    min=0,
    max=3,
    step=0.01,
    description='alpha:',
    disabled = False ,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to parallax alpha-slider

alpha = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((alpha0, 'value'), (alpha, 'value'))

#beta-slider

beta0 = widgets.FloatSlider(
    value=2,
    min=0,
    max=3.00,
    step=0.01,
    description='beta:',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to slider

beta = widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((beta0, 'value'), (beta, 'value'))

#mu_ra

mu_ra0 = widgets.FloatSlider(
    value=0,
    min=-60,
    max=60,
    step=0.1,
    description='propm ra [mas/yr]: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to slider

mu_ra= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((mu_ra0, 'value'), (mu_ra, 'value'))

#mu_dec

mu_dec0 = widgets.FloatSlider(
    value=0,
    min=-60,
    max=60,
    step=0.1,
    description='propm dec [mas/yr]: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to slider

mu_dec= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((mu_dec0, 'value'), (mu_dec, 'value'))

#sd_mu_ra

sd_mu_ra0 = widgets.FloatSlider(
    value=1,
    min=-60,
    max=60,
    step=0.1,
    description='propm ra error[mas/yr]: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)

#textfield connected to slider

sd_mu_ra= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((sd_mu_ra0, 'value'), (sd_mu_ra, 'value'))


#sd_mu_dec

sd_mu_dec0 = widgets.FloatSlider(
    value=1,
    min=-60,
    max=60,
    step=0.1,
    description='propm dec error[mas/yr]: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)
sd_mu_dec= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((sd_mu_dec0, 'value'), (sd_mu_dec, 'value'))



#corr_w_mu_ra

corr_w_mu_ra0 = widgets.FloatSlider(
    value=0,
    min=-1,
    max=1,
    step=0.1,
    description='parallax - propm ra correlation: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)
corr_w_mu_ra= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((corr_w_mu_ra0, 'value'), (corr_w_mu_ra, 'value'))

#corr_w_mu_dec
corr_w_mu_dec0 = widgets.FloatSlider(
    value=0,
    min=-1,
    max=1,
    step=0.1,
    description='parallax - propm dec correlation: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)
corr_w_mu_dec= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((corr_w_mu_dec0, 'value'), (corr_w_mu_dec, 'value'))

#corr_mu_ra_dec

corr_mu_ra_dec0 = widgets.FloatSlider(
    value=0,
    min=-1,
    max=1,
    step=0.1,
    description='propm ra - propm dec correlation: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)
corr_mu_ra_dec= widgets.FloatText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((corr_mu_ra_dec0, 'value'), (corr_mu_ra_dec, 'value'))

#healpix

healpix0 = widgets.IntSlider(
    value=6200,
    min=0,
    max=12288,
    step=1,
    description='healpix: ',
    disabled = False,
    continuous_update=False,
    orientation='horizontal',
    readout=False,
    #readout_format='.2f',
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)
healpix= widgets.IntText(layout=widgets.Layout(width='25%'),flex='2')
widgets.jslink((healpix0, 'value'), (healpix, 'value'))

healpix_description=widgets.Label(value='only used for velocity prior; not for geometric distance prior')



## walker initialisation
#
#walker_init = widgets.Text(
#    value='0.31, 0.35, 0.36, 0.37',
#    description='walker initialisation [kpc]',
#    disabled=False,
#    style = {'description_width': 'initial'},
#    layout=widgets.Layout(flex='8')
#)

## plotting range:
#
#plotting_range = widgets.FloatRangeSlider(
#    value=[0, 5],
#    min=0,
#    max=12,
#    step=0.01,
#    description='Plotting range:',
#    disabled=False,
#    continuous_update=False,
#    orientation='horizontal',
#    style = {'description_width': 'initial'},
#    layout=widgets.Layout(flex='8')
#
#)

#create pdf option

create_pdf = widgets.Checkbox(
    value=False,
    description='Create .pdf file containing plots for all sources',
    disabled=False,
    style = {'description_width': 'initial'},
    layout=widgets.Layout(flex='8')
)


# display all widgets

display(sampler)
display(data)

display(widgets.HBox([source_id, description_text]))
display(Markdown('(Mode "single source, source_id" must be selected to use this.)'))
display(filename_in)
display(filename_out)
display(create_pdf)
display(seed)
display(Markdown(" "))
display(Markdown("**General sampling settings:** "))
display(Markdown("(are used for each possible mode) "))
hbox_layout = widgets.Layout(display='flex', flex_flow='row', justify_content='space-between', width='100%')

display(widgets.HBox([Nsamp0,Nsamp],layout=hbox_layout))
display(widgets.HBox([Nburnin0,Nburnin],layout=hbox_layout))
display(widgets.HBox([thinfac0,thinfac],layout=hbox_layout))
display(widgets.HBox([n0,n,widgets.Label(value="(set to 1 if unsure what this is)")],layout=hbox_layout))
#display(Markdown("(set to 1 if unsure what this is)"))
display(Markdown(" "))
display(Markdown("**Sampling settings only required for emcee:** "))
display(widgets.HBox([Nwalker0,Nwalker],layout=hbox_layout))
display(Markdown ('Initialisation for walkers: random values between (1-a)*starting value to (1+a)*starting value:'))
display(widgets.HBox([a0,a],layout=hbox_layout))
display(Markdown(" "))
display(Markdown("**Sampling settings only reqired when using own data for a single source:**"))
display(widgets.HBox([rInit0,rInit],layout=hbox_layout))
display(Markdown ('(Initialisation when using source_ids: mode of EDSD prior)'))
display(widgets.HBox([rStep0,rStep],layout=hbox_layout))
#display(Markdown("**Sampling settings only required for metropolis when processing a single source:**"))
#display(Markdown("(when processing multiple sources, the initialisation will be 1/abs(parallax) and the stepsize will be 0.75* initialisation* min(1/3,abs(parallax_error/parallax)))"))

#display(Markdown(" "))
#display(Markdown("**Sampling settings only required for emcee when processing a single source:** "))

# display input data 
display(Markdown(" "))
display(Markdown("**Input data when using own input data for a single source:**"))
display(Markdown("(Default alpha and beta values are for EDSD prior)"))
display(widgets.HBox([parallax0,parallax],layout=hbox_layout))
display(widgets.HBox([parallax_error0,parallax_error],layout=hbox_layout))
display(widgets.HBox([rlen0,rlen],layout=hbox_layout))
display(widgets.HBox([alpha0,alpha],layout=hbox_layout))
display(widgets.HBox([beta0,beta],layout=hbox_layout))
display(widgets.HBox([mu_ra0,mu_ra],layout=hbox_layout))
display(widgets.HBox([mu_dec0,mu_dec],layout=hbox_layout))
display(widgets.HBox([sd_mu_ra0,sd_mu_ra],layout=hbox_layout))
display(widgets.HBox([sd_mu_dec0,sd_mu_dec],layout=hbox_layout))
display(widgets.HBox([corr_w_mu_ra0,corr_w_mu_ra],layout=hbox_layout))
display(widgets.HBox([corr_w_mu_dec0,corr_w_mu_dec],layout=hbox_layout))
display(widgets.HBox([corr_mu_ra_dec0,corr_mu_ra_dec],layout=hbox_layout))
display(widgets.HBox([healpix0,healpix,healpix_description],layout=hbox_layout))
display(Markdown(" "))

#display(Markdown("**Plotting range when processing a single source:**"))
#display(Markdown('(when processing multiple sources, the range will be 0.2*min(kinegeo samples)-1.2*max(kinegeo samples))'))
#display(widgets.VBox([plotting_range],layout=hbox_layout))

out = widgets.Output()

# submit button

submit_button = widgets.Button(description='start')
display(submit_button)

# submit function: when clicking submit, the simulation is run for all settings that have been made

def submit(button):
    out.clear_output()
    with out:
        interactive_velocity_function(data = data.value, sampler=sampler.value, seed=seed.value, source_id = source_id.value, 
                                      parallax=parallax.value, parallax_error=parallax_error.value, 
                                      alpha=alpha.value, beta=beta.value, rlen=rlen.value, 
                                      mu_ra=mu_ra.value, mu_dec=mu_dec.value, sd_mu_ra=sd_mu_ra.value, sd_mu_dec=sd_mu_dec.value, 
                                      corr_w_mu_ra= corr_w_mu_ra.value, corr_w_mu_dec=corr_w_mu_dec.value,
                                      corr_mu_ra_dec=corr_mu_ra_dec.value, 
                                      healpix=int(healpix.value),
                                      Nwalker=Nwalker.value, a=a.value, 
                                      rInit=rInit.value, rStep=rStep.value, Nsamp=Nsamp.value, 
                                      thinfac=thinfac.value, Nburnin=Nburnin.value, n=n.value, 
                                      filename_in = filename_in.value,filename_out = filename_out.value, create_pdf = create_pdf.value,
                                      rows_prior_summary=rows_prior_summary, Nmax=Nmax, probs = probs)
        
# tie submit button to a function
submit_button.on_click(submit)

out 
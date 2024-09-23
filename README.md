# Interactive-Kinegeometric-Distance-and-Velocity-Estimation

The Jupyter Notebook "Interactive Kinegeometric Distance and Velocity Estimation.ipynb" can be used to infer combined distances and velocities using a kinegeometric distance and velocity posterior. There are two different sampling methods to choose from: the metropolis algorithm and the python emcee package.

To start the interactive distance and velocity estimation, run all cells first. Then choose your settings and press 'start' to begin with the distance estimation. At first you choose the sampler you want to use, there is either the metropolis or the emcee option. Then you choose the mode. At first there is the option to process a single source using own input data, which can be specified using the sliders and textfields below. Then there is the option to either insert a Gaia DR3 source_id or the name of a specific star into the field below. For this you have to select the mode 'single source, source_id'. All the data for this source is then retrieved from Gaia DR3. 

The last option does the same, but for a .csv file of multiple source_ids specified in the 'Input .csv-file' field. The results are then saved to three .csv files in the 'results' folder, which can be named in the 'Output .pdf/.csv-file' field. There is the  '[custom name]_samples_geo.csv' file, which contains the geometric samples for each source. Then there is the '[custom name]_samples_kinegeo.csv' file, which contains the kinegeometric distance samples and finally, there is the '[custom_name]_summary_statistics.csv' file containing all summary statistics. Additionally, a .pdf file containing distance and velocity plots of all sources can be created when ticking the checkbox 'Create .pdf file containing plots for all sources'. 

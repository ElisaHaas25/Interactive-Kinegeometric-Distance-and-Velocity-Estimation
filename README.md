# Interactive-Kinegeometric-Distance-and-Velocity-Estimation

The Jupyter Notebook "Interactive Kinegeometric Distance and Velocity Estimation.ipynb" can be used to infer combined distances and velocities using a kinegeometric distance and velocity posterior. There are two different sampling methods to choose from: the metropolis algorithm and the python emcee package.

Open the notebook ‘Interactive Kinegeometric Distance Estimation.ipynb. To start the interactive distance and velocity estimation, run the cell. Then choose your settings and press 'start' to begin with the distance estimation. At first you choose the sampler you want to use, there is either the metropolis or the emcee option. Then you choose the mode. At first there is the option to process a single source using own input data, which can be specified using the sliders and textfields below. Then there is the option to either insert a Gaia DR3 source_id or the name of a specific star into the field below. For this you have to select the mode 'single source, source_id'. All the data for this source is then retrieved from Gaia. 
The last option does the same, but for a .csv file of multiple source_ids specified in the 'Input .csv-file' field. Here you can enter the name of your file without the '.csv' suffix. The file has to be in the 'data' folder. 

For all Modes, the results are saved to four .csv files in the 'results' folder, which can be named in the 'Output .pdf/.csv-file' field. There is the  '[custom name]_distanceSamplesGeo.csv' file, which contains the geometric samples for each source. There is the '[custom name]_distanceSamplesKinegeo.csv' file, which contains the kinegeometric distance samples. Then there is the '[custom name]_velocitySamples.csv' file, which contains the Ra and Dec velocity samples and finally, there is the '[custom_name]_summaryStatistics.csv' file containing all summary statistics. Additionally, a .pdf file containing distance and velocity plots of all sources can be created when ticking the checkbox 'Create .pdf file containing plots for all sources'. 

The estimation of kinegeometric distances is described in the paper
“Estimating distances from parallaxes. VI. A method for inferring distances and transverse velocities from parallaxes and proper motions demonstrated on Gaia Data Release 3” by Coryn Bailer-Jones in 2023 (Astronomical Journal, 166, 269)
https://doi.org/10.3847/1538-3881/ad08bb
For more details see: https://www.mpia.de/homes/calj/gedr3_distances_velocities.html
The method implemented in the current notebook is slightly different, and is described here:
https://www.mpia.de/homes/calj/gedr3_distances_velocities/distance_velocity_sampling_1d.pdf


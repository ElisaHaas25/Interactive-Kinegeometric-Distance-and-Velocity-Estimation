# load all required packages

#general packages 
import numpy as np
import math
import re
import matplotlib.pyplot as plt
import csv

#scipy 
from scipy import integrate
from scipy.stats import norm
from scipy.special import gamma, factorial, loggamma
from scipy.stats import multivariate_normal as mvn

#rpy2 
import rpy2
from rpy2.robjects.packages import importr, data
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
from rpy2.robjects import r
numpy2ri.activate() #activate conversion from R-code to numpy arrays

# gaia, zeropoint correction
from astroquery.gaia import Gaia
from zero_point import zpt # need to install gaiadr3-zeropoint
zpt.load_tables()
from astroquery.simbad import Simbad

#samplers 
import emcee #https://emcee.readthedocs.io/en/stable/
from metropolis import metrop
#from functions import mode_post3 

#for interactive display
import ipywidgets as widgets
from IPython.display import display, Markdown
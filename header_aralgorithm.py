### NETCDF/datetime stuff
from netCDF4 import Dataset,num2date,date2num
import datetime
from datetime import timedelta, datetime, date


### Matplotlib stuff
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.transforms as mtrans
import matplotlib.patheffects as pe
from matplotlib import cm


### Data stuff
import numpy as np
import numpy.ma as ma
from scipy import interpolate
import sys, os, os.path
import pandas as pd


### Mapping stuff
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean


#### AR detection algorithm dependencies
import warnings
warnings.filterwarnings('ignore')
# Import modules and packages
# import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import csv
# from datetime import date
from dateutil.relativedelta import relativedelta
import geopy.distance
import iris
import iris.coords
import iris.quickplot as qplt
# import matplotlib
# import matplotlib.pyplot as plt
# import numpy as np
import operator
# import os
# import os.path
from scipy import ndimage
from shapely.geometry import Point
# import skimage
# from skimage.segmentation import find_boundaries
import xarray as xr
import xesmf 


### Utility stuff
import timeit
import gc
import threading as tr
import json
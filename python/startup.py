import sys,os,glob
import subprocess

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LogNorm,Normalize,SymLogNorm
import cmocean.cm as cmo
import cmasher as cma

import pandas as pd
import xarray as xr

import astropy.constants as ac
import astropy.units as au

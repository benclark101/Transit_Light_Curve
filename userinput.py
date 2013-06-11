import pyfits
import numpy
import csv
import os
import math
import matplotlib.pyplot as plt
from scipy import special as sp
import random
import scipy.optimize as optimize
import scipy.integrate as intg

global fitsLoc
global OutFolder 
global IntenOutNm 
global ErrOutNm 
global rMax 
global pi
global a
global rsun

#--------------- VARIABLES ------------------
fitsLoc = 'C:/Users/ben/Desktop/Project - 4th Year/20120101/'#location of fits files
#Ensure that only fits files are in the folder given above!
OutFolder = 'Output/' #Folder to Output data
#--------------Constants-----------------
rMax = 14 #Aperture radius size
pi = 3.14159265358979323846264338327 #pi
a = 0
rsun = 6.955E8 # In m
IntenOutNm = 'data.txt' #DO NOT CHANGE
ErrOutNm = 'error.txt' #DO NOT CHANGE
#---------Host Star Properties-------------
Mstar = 1.35 * 2E30 #in KG
ErMstar = 0.14 * 2E30 #in KG
Rstar = 1.57 * 6.955E8 # In m
ErRstar = 0.07 * 6.955E8 # In m


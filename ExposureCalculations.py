"""
This module is used to calculate exposure times according to recipe from Burt et al. (2015).
The basic procedure is as follows:
1. For a specific star (V mag and B-V color) a precision is requested.
2. That precision and the color in turn generate a required median number of counts in the I2 region (getI2_M or getI2_K), the cutoff betweeen K and M is B-V > 1.2 for M
3. Given the color, the magnitude, the seeing (in pixels on the guider) and the elevation, the rate of photons landing in the I2 region of spectrometer can be calculated.
4. The ratio of line 2 to line 3 is the exposure time.
5. A certain number of I2 counts corresponds to a certain number of exposure meter counts, which has a color correction. This is the second possible way to control the epxosure time.

A number of modules rely on the dictionary slit_size which is defined at the top of the module. This has the slit dimensions a dictionary of tuples. The dimensions are in order of (width, length) in pixels. So the W decker, which is 1" by 3" is (9.26, 27.78).

"""

import numpy as np
from x_gaussslit import *

# slit_size is a dictionary that defines the apertures of the APF in pixels. 

slit_size = {'M': (9.26, 74.07), 'W': (9.26, 27.78), 'N': (4.63, 74.07), 'B': (18.52, 74.07), 'S': (6.94, 74.07), 'P': (1, 1), 'O' : (74.07,74.07), 'T': (18.52, 27.78),  'L': (18.52, 111.11), 'K': (9.26, 111.11)}
    
def getI2_K(unc):
    """
    i2counts = getI2_K(unc)
    unc - uncertainty in m/s

    i2counts - median number of counts in the I2 region of the spectrum.

    This relation is for K and G stars, B-V < 1.2.
    """
    A = 4.47
    B = -1.58
    return 10 ** ( A + B*np.log10(unc))
    
def getI2_M(unc):
    """
    i2counts = getI2_M(unc)
    unc - uncertainty in m/s

    i2counts - median number of counts in the I2 region of the spectrum.

    This relation is for K and G stars, B-V > 1.2.
    """


    A = 4.14
    B = -1.73
    return 10 ** ( A + B*np.log10(unc))
    
def getEXPMeter(i2, bv):
    """
    expcounts = getEXPMeter(i2counts, bmv)

    expcounts - number of counts in the exposure meter when a star of the input color has the input number of counts in the I2 region

    i2counts - number of counts desired in the spectrum
    bv - color of the star (B-V, Johnson filters, Vega mags.)
    """
    delta = 4.52
    epsilon = -0.196
    squiggle = 0.262
    x = delta + epsilon*bv + squiggle*bv**2
    x += 0.05
    return i2 * 10**x

def getEXPMeter_Rate(v, bv, el, seeing, decker="W"):
    """
    expcountrate = getEXPMeter_Rate(vmag, bmv, elevation, seeing, decker="W")

    expcountrate - rate of photons landing on the exposure meter
    
    vmag - V magnitude of star (Johnson filters, Vega mag)
    bv - color of the star (B-V, Johnson filters, Vega mags.)
    el - elevation in degrees
    seeing - seeing in pixels on the guider
    decker - an ascii letter, must match value in dictionary defined at top of module
    """
    alpha = -0.908
    beta = 0.0852
    Const = -21.8
    if seeing == 0:
#        apflog( "Warning: AVG_FWHM seems to be 0. Using 15 instead.",level="Warn")
        seeing = np.array(15)
    # seeing  = 13.99
    light = x_gaussslit(slit_size[decker][0]/seeing, slit_size[decker][1]/seeing, 0, 0)
    # light = 0.442272

    if light > 0:    
        VC = v - 2.5*np.log10(light)
        x = (-1/2.5) * (VC + alpha*bv + beta*(1/np.cos(np.radians(90-el))) + Const)
    else:
        x = -10
    return (10 ** x)

def getSpec_Rate( v, bv, el, seeing, decker="W"):
    """
    rate = getSpec_Rate(vmag, bmv, el, seeing, decker="W")

    rate - arrival rate in photons per second in the I2 region of the spectrometer.
    
    vmag - V magnitude of star (Johnson filters, Vega mag)
    bv - color of the star (B-V, Johnson filters, Vega mags.)
    el - elevation in degrees
    seeing - seeing in pixels on the guider
    decker - an ascii letter, must match value in dictionary defined at top of module
    """
    alpha = -0.0311
    beta = 0.158
    Const = -11.2
    if seeing <= 0:
        seeing = np.array(15)
    # seeing  = 13.99
    light = x_gaussslit(slit_size[decker][0]/seeing, slit_size[decker][1]/seeing, 0, 0)
    # light = 0.442272
    

    if light > 0:
        VC = v - 2.5*np.log10(light)
        x = (-1/2.5) * (VC + alpha*bv + beta*(1/np.cos(np.radians(90-el))) + Const)
    else:
        x = -10
    cnt_rate = 10**x
    return cnt_rate

def getEXPTime(cnts, v, bv, el, seeing, decker="W"):
    """ time = getEXPTime(cnts, v, bv, el, seeing, decker="W")

    time - time in seconds
    
    vmag - V magnitude of star (Johnson filters, Vega mag)
    bmv - color of the star (B-V, Johnson filters, Vega mags.)
    el - elevation in degrees
    seeing - seeing in pixels on the guider
    decker - an ascii letter, must match value in dictionary defined at top of module

    """
    time = 0
    cnt_rate = getSpec_Rate(v, bv, el, seeing, decker="W")
    fin_cnt_rate = np.where(cnt_rate > 0,cnt_rate,1.0e-5)
    time = cnts/fin_cnt_rate
    return time

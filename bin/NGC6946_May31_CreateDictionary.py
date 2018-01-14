__author__ = 'Locke Patton'
# Run 'ulimit -n 512' before executing

#CHANNELS - EMISSION LINES

# blue channel
# [OII] 3727*    *E
# [OIII] 4363 (weak)
# He II 4686 (weak)
# Hbeta 4861*    *E
# [OIII] 4959
# [OIII] 5007*    *E

# red channel
# [OI] 6300 (weak)
# [NII] 6548
# Halpha 6563*
# [NII] 6584 (weak)
# [SII] 6717
# [SII] 6731

from SlitSpectra5 import SlitSpectra as SS
from SlitSpectra5 import findRADEC
import warnings
import os
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.pyplot import cm
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import astropy.coordinates as coord
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle,Rectangle
from astropy.visualization import astropy_mpl_style
from astropy.table import Table, hstack
from PyAstronomy import pyasl
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import mode

# Verbose
verbose = False
verbose_absolute = True

# Directories
#line profile inputs
#location of slot.log files
dir1_ = r'/Users/lockepatton/Desktop/Research/Levesque/NGC6946/data2.0/'
dir2_ = r'apall4_x2.0_fluxcalibration/4.2calibrate/x_dFT_red/'
dir_lineprofile_ = dir1_ + dir2_

# y dir
# '/Users/lockepatton/Desktop/Research/Levesque/NGC6946/data2.0/' + \
# 'apall4_y_fluxcalibration/apall4_y_after_4.2fluxcalibration/'
# x dir
# '/Users/lockepatton/Desktop/Research/Levesque/NGC6946/data2.0/' + \
# 'apall4_x2.0/2files_dFT_withbackground/'

#defining Images
dTImages_r = ['dT4.1_02hh_NS.0052r',
            'dT4.1_04et_NS.0047r',
            'dT4.1_08S_NS.0051r',
            'dT4.1_17A_NS.0050r',
            'dT4.1_39C_NS.0053r',
            'dT4.1_48B_NS.0046r',
            'dT4.1_68D_NS.0045r',
            'dT4.1_69P_NS.0049r',
            'dT4.1_80K_NS.0048r']
dFTImages_r = ['dFT4.1_02hh_NS.0052r',
            'dFT4.1_04et_NS.0047r',
            'dFT4.1_08S_NS.0051r',
            'dFT4.1_17A_NS.0050r',
            'dFT4.1_39C_NS.0053r',
            'dFT4.1_48B_NS.0046r',
            'dFT4.1_68D_NS.0045r',
            'dFT4.1_69P_NS.0049r',
            'dFT4.1_80K_NS.0048r']
cdFTImages_r = ['cdFT4.1_02hh_NS.0052r',
            'cdFT4.1_04et_NS.0047r',
            'cdFT4.1_08S_NS.0051r',
            'cdFT4.1_17A_NS.0050r',
            'cdFT4.1_39C_NS.0053r',
            'cdFT4.1_48B_NS.0046r',
            'cdFT4.1_68D_NS.0045r',
            'cdFT4.1_69P_NS.0049r',
            'cdFT4.1_80K_NS.0048r']
dTImages_b = ['dT4.1_02hh_NS.0052b',
            'dT4.1_04et_NS.0047b',
            'dT4.1_08S_NS.0051b',
            'dT4.1_17A_NS.0050b',
            'dT4.1_39C_NS.0053b',
            'dT4.1_48B_NS.0046b',
            'dT4.1_68D_NS.0045b',
            'dT4.1_69P_NS.0049b',
            'dT4.1_80K_NS.0048b']
dFTImages_b = ['dFT4.1_02hh_NS.0052b',
            'dFT4.1_04et_NS.0047b',
            'dFT4.1_08S_NS.0051b',
            'dFT4.1_17A_NS.0050b',
            'dFT4.1_39C_NS.0053b',
            'dFT4.1_48B_NS.0046b',
            'dFT4.1_68D_NS.0045b',
            'dFT4.1_69P_NS.0049b',
            'dFT4.1_80K_NS.0048b']
cdFTImages_b = ['cdFT4.1_02hh_NS.0052b',
            'cdFT4.1_04et_NS.0047b',
            'cdFT4.1_08S_NS.0051b',
            'cdFT4.1_17A_NS.0050b',
            'cdFT4.1_39C_NS.0053b',
            'cdFT4.1_48B_NS.0046b',
            'cdFT4.1_68D_NS.0045b',
            'cdFT4.1_69P_NS.0049b',
            'cdFT4.1_80K_NS.0048b']

def trace_center_r_or_b(channel, trace_center_r, trace_center_b):
    if channel == 'b':
        trace_center = trace_center_b
    if channel == 'r':
        trace_center = trace_center_r
    return trace_center

def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

#function that creates a SlitSpectra Object in X + Y directions
#Building Y direction
def CreateSS_y(image,plotSpec=False,verbose=False):
    channel = image[-1]
    dir_ = '/Users/lockepatton/Desktop/Research/Levesque/NGC6946/data2.0/apall4_y_fluxcalibration/apall4_y_after_4.2fluxcalibration/'
    full_region = [1,2098,1,1078]
    center_width = 10
    trace_center = 0 #initially set before know value.

    trim_region_r=[15,1920,99,819]
    full_lowest_r=83.5
    full_highest_r=834.5

    trim_region_b=[210,1990,210,888]
    full_lowest_b=195
    full_highest_b=903.5

    if channel == 'r':
        trace_center = 65.29
        trim_region = trim_region_r
        full_lowest = full_lowest_r
        full_highest = full_highest_r
        w1 =  5428.431 #starting wavelength
        w2 =  9802.673 #ending wavelength
        dw =  2.29619 #wavelength intervel per pixel
        nw =  1906 #number of output pixels

    if channel == 'b':
        trace_center = 1458.02
        trim_region = trim_region_b
        full_lowest = full_lowest_b
        full_highest = full_highest_b
        w1 =  2926.267 #starting wavelength
        w2 =  6165.38 #ending wavelength
        dw =  1.819726 #wavelength intervel per pixel
        nw =  1781 #number of output pixels

    trim_base_vs_transfer_fit = [0.94335225214255913, 0.81555279062263253]
    c1,c2 = trim_base_vs_transfer_fit

    T = SS(image=image,
           channel=channel,
           path=dir_,
           full_region=full_region,
           trim_region=trim_region,
           direction='y',
           full_lowest=full_lowest,
           full_highest=full_highest,
           trace_center=trace_center,
           SN_width=center_width,
           w1=w1,w2=w2,dw=dw)

    wavelength_widths = [6,6,6,6,6,6]

    if channel == 'b':
        wavelengths = [3727, 4363, 4686, 4861, 4959, 5007]   #"blue" wavelengths
        wavelength_names = [r'[OII]3727*', r'[OIII]4363weak',
                            r'HeII4686weak', r'Hbeta4861*',
                            r'[OIII]4959', r'[OIII]5007*'] #first half
        wavelength_widths = [6,6,6,6,6,6]

    if channel == 'r':
        wavelengths = [6300, 6548, 6563, 6584, 6717, 6731]     #red wavelengths
        wavelength_names = [r'[OI]6300weak', r'[NII]6548',
                            r'Halpha6563*', r'[NII]6584weak',
                            r'[SII]6717', r'[SII]6731']
        wavelength_widths = [6,6,6,6,6,6]

    pixel_backgrounds = []

    it_ = 0
    while it_ < len(wavelengths):
        pixel_backgrounds.append([-7,-5,5,7])
        it_ = it_ + 1

    T.calculateWavelengths(wavelength_names=wavelength_names,
                           wavelengths=wavelengths, pixel_widths=wavelength_widths,
                           pixel_backgrounds=pixel_backgrounds)

    if plotSpec:
        T.plotSpectra(save_fig=False, verbose=False)
        # T.plotImageWavelength(save_fig=False,verbose=False)

    T.buildSpectra(verbose=verbose)
    return T

#Building X objects
def CreateSS_x(image,plotSpec=False,verbose=False):
    channel = image[-1]
    dir_ = '/Users/lockepatton/Desktop/Research/Levesque/NGC6946/data2.0/apall4_x2.0/2files_dFT_withbackground/'
    full_region = [1,2098,1,1078]
    center_width = 10
    trace_center = 0       #initially set before know value.

    trim_region_r=[15,1920,99,819]
    full_lowest_r=83.5
    full_highest_r=834.5

    trim_region_b=[210,1990,210,888]
    full_lowest_b=195
    full_highest_b=903.5

    if channel == 'r':
        trim_region = trim_region_r
        full_lowest = full_lowest_r
        full_highest = full_highest_r
        w1 =  5428.431 #starting wavelength
        w2 =  9802.673 #ending wavelength
        dw =  2.29619 #wavelength intervel per pixel
        nw =  1906 #number of output pixels

    if channel == 'b':
        trim_region = trim_region_b
        full_lowest = full_lowest_b
        full_highest = full_highest_b
        w1 =  2926.267 #starting wavelength
        w2 =  6165.38 #ending wavelength
        dw =  1.819726 #wavelength intervel per pixel
        nw =  1781 #number of output pixels

    trim_base_vs_transfer_fit = [0.94335225214255913, 0.81555279062263253]
    c1,c2 = trim_base_vs_transfer_fit

    if image[-17:-1] == '4.1_02hh_NS.0052':
        trace_center_r = 489
        trace_center_b = 139.63
    if image[-17:-1] == '4.1_04et_NS.0047':
        trace_center_r = 159.12 #157.775
        trace_center_b = 148
    if image[-16:-1] == '4.1_08S_NS.0051':
        trace_center_r = 712.207
        trace_center_b = 670
    if image[-16:-1] == '4.1_17A_NS.0050':
        trace_center_r = 313.38
        trace_center_b = 294
    if image[-16:-1] == '4.1_39C_NS.0053':
        trace_center_r = 652.623
        trace_center_b = 612
    if image[-16:-1] == '4.1_48B_NS.0046':
        trace_center_r = 501.5273 #484
        trace_center_b = 467.31
    if image[-16:-1] == '4.1_68D_NS.0045':
        trace_center_r = 288.14  #287.669
        trace_center_b = 238
    if image[-16:-1] == '4.1_69P_NS.0049':
        trace_center_r = 427.26 #426.928
        trace_center_b = 402
    if image[-16:-1] == '4.1_80K_NS.0048':
        trace_center_r = 487.76 #488.154
        trace_center_b = 230

    trace_center = trace_center_r_or_b(channel, trace_center_r, trace_center_b)

    T = SS(image=image,
           channel=channel,
           path=dir_,
           full_region=full_region,
           trim_region=trim_region,
           direction='x',
           full_lowest=full_lowest,
           full_highest=full_highest,
           trace_center=trace_center,
           SN_width=center_width,
           w1=w1,w2=w2,dw=dw)

    T.findcenterwidthBackground(verbose=verbose)
    T.allcalc(verbose=verbose)

    if plotSpec:
    #     T.plotSpectra(x=[6000,7000],save_fig=True,verbose=False)
    #     T.plotSpectra(x=[4750,5750],save_fig=True,verbose=False)

        T.plotSpectra(verbose=False)

    #     T.plotImage(save_fig=False,vmin=0, vmax=250, verbose=True, blackout=True)
    # #     T.plotImage(save_fig=False, x=[250,2500], verbose=False, blackout=False)
    #     if T.channel == 'b':
    #         print 'plotImage', T.channel
    #         T.plotImage(x=[1000,1200],save_fig=False,verbose=False)  #b Signal in b
    #     if T.channel == 'r':
    #         print 'plotImage', T.channel
    #         T.plotImage(x=[330,700],save_fig=False,verbose=False)  #r Signal in r - Halpha Group

    return T

#Defining N2, log(O/H)+12 polynomial fit and log(O/H)+12 linear fit functions
func_N2 = lambda NII,Halpha : np.log10(NII/Halpha)
func_logOH = lambda N2 : 9.37 + 2.03*N2 + 1.26*N2**2 + 0.32*N2**3
func_logOH_linear = lambda N2 : 8.90 + 0.57* N2

#CREATING DICTIONARY TO SAVE
Met_N2 = {}
Met_N2['image'] = []
Met_N2['x_pix'] = []
Met_N2['AllSpec'] = []
Met_N2['Halpha6563*'] = []
Met_N2['[NII]6584weak'] = []
Met_N2['Halpha_bkgdsub'] = []
Met_N2['NII_bkgdsub'] = []
Met_N2['Halpha_sigmas']=[]
Met_N2['NII_sigmas'] =[]
Met_N2['N2'] = []
Met_N2['logOH_linear']  = []
Met_N2['logOH'] = []
Met_N2['HighSignalRegions'] = []
Met_N2['StarRegions'] = []
Met_N2['Backgrounds'] = {}
Met_N2['Backgrounds']['mean'] = {}
Met_N2['Backgrounds']['mode'] = {}
Met_N2['Backgrounds']['sigma_bkgdsub'] = {}
Met_N2['Backgrounds']['sigma_bkgdsub']['Ha'] = []
Met_N2['Backgrounds']['sigma_bkgdsub']['NII'] = []
Met_N2['Backgrounds']['regions_bkgdsub'] = {}
Met_N2['Backgrounds']['regions_bkgdsub']['Ha'] = []
Met_N2['Backgrounds']['regions_bkgdsub']['NII'] = []
Met_N2['Backgrounds']['background_fits_all'] = []
Met_N2['AllSpec_Background_Mode_Dict'] = {}
Met_N2['RA'] = []
Met_N2['DEC'] = []
Met_N2['RA_array'] = []
Met_N2['DEC_array'] = []
Met_N2['RAh'] = []
Met_N2['RAm'] = []
Met_N2['RAs'] = []
Met_N2['DECd'] = []
Met_N2['DECm'] = []
Met_N2['DECs'] = []
Met_N2['NWSE'] = []


for dT_image_r, cdFTImage_r, dFT_image_r in zip(dTImages_r, cdFTImages_r, dFTImages_r):

    image = dFT_image_r
    S = CreateSS_x(image)
    # S.buildSpectra()  # don't believe it's necissary to build spectra in x direction

    image = cdFTImage_r
    T = CreateSS_y(image,plotSpec=False)
    T.buildSpectra()

    returnsAllSpec = T.returnAllSpectra()
    [x_pix, SpectraAll] = returnsAllSpec

    # returns = T.plotWavelengthSpectra('[NII]6584weak','Halpha6563*',x=None,minSubtract=False,
    #                         verbose=False,returned=True,save_fig=True)
    #
    # [x_pix, [y_NII6584_min_sub, y_NII6584_true, wavelength_NII6584], \
    #         [y_Halpha6563_min_sub, y_Halpha6563_true, wavelength_Halpha6563]] = returns

    Met_N2['image'].append(image)
    Met_N2['AllSpec'].append(SpectraAll)
    Met_N2['x_pix'].append(x_pix)
    Met_N2['Halpha6563*'].append(SpectraAll['Halpha6563*'])
    Met_N2['[NII]6584weak'].append(SpectraAll['[NII]6584weak'])
    del T


#High Signal Regions Within Each Image
dictionary_no_regions = {
    'cdFT4.1_02hh_NS.0052r' : [[410, 520], [590, 610]], #1 cutting only region w\ 2-prong spikes, and bad pixel in NII
    'cdFT4.1_04et_NS.0047r' : [[0, 90], [140, 180], [230, 260]], #2 three very strong regions of high spikes
    'cdFT4.1_08S_NS.0051r' : [[40, 90], [180, 300], [580, 600], [690, 721]], #3 three major regions with signal spikes
    'cdFT4.1_17A_NS.0050r' : [[40, 80], [290, 560], [590, 721]], #4 tight cutting, very little lest
    'cdFT4.1_39C_NS.0053r' : [[300, 460], [530, 721]], #5 tight cutting. only two 'background' regions left
    'cdFT4.1_48B_NS.0046r' : [[420, 550], [660, 700]], #6 took out 2 major 'high' regions that hit ~.25 in Halpha
    'cdFT4.1_68D_NS.0045r' : [[180, 300], [370, 400], [430, 450], [500, 620]], #7 lots of small-ish signal spikes removed
    'cdFT4.1_69P_NS.0049r' : [[0, 130], [180, 270], [410, 460], [540, 590]], #8 3 major spikes and spikey reigon.. left inverse spike
    'cdFT4.1_80K_NS.0048r' : [[0, 130], [230, 260], [470, 500]] #9 took out higher edge region in low xpix, and two high spikes
}

#Stars + Bad Pixels
dictionary_star_regions = {
    'cdFT4.1_02hh_NS.0052r' : [[595,610]],  #bad pixel region
    'cdFT4.1_04et_NS.0047r' : [[150,170],[240,253]], #def 2 stars
    'cdFT4.1_08S_NS.0051r' : [[585,590]], #bad pixel region
    'cdFT4.1_17A_NS.0050r' : [],
    'cdFT4.1_39C_NS.0053r' : [],
    'cdFT4.1_48B_NS.0046r' : [],
    'cdFT4.1_68D_NS.0045r' : [[283,293]], #might be a star - check - NII ~ Halpha
    'cdFT4.1_69P_NS.0049r' : [[421,433]], #dim but star
    'cdFT4.1_80K_NS.0048r' : [[480,493]], #star
}

Met_N2['no_regions'] = dictionary_no_regions
Met_N2['star_regions'] = dictionary_star_regions


# All running background regions - new code (keeping old code for backward compatability)
wavelength_profiles = [r'[OI]6300weak', r'[NII]6548',
                       r'Halpha6563*', r'[NII]6584weak',
                       r'[SII]6717', r'[SII]6731']
for specname in wavelength_profiles:
    MassiveBkgdSpec = []
    for image, AllSpec in zip(Met_N2['image'],Met_N2['AllSpec']):
        spectra = AllSpec[specname].copy()
        for no_region in dictionary_no_regions[image]:
            if len(no_region) > 0:
                no_region_low, no_region_hi = no_region
                spectra[no_region_low:no_region_hi] = np.nan
        MassiveBkgdSpec.append(spectra)

    MassiveBkgdSpec = np.array(MassiveBkgdSpec)
    Mode = mode(MassiveBkgdSpec)[0][0]
    Met_N2['AllSpec_Background_Mode_Dict'][specname] = Mode


# Building complete 'background' low signal regions
T_xpix = []
T_Ha = []
T_NII = []

for i, [image,x_pix, Ha, NII] in enumerate(zip(Met_N2['image'],Met_N2['x_pix'],Met_N2['Halpha6563*'],Met_N2['[NII]6584weak'])):
    Ha_backgroundonly = Ha.copy()
    NII_backgroundonly = NII.copy()
    for no_region in dictionary_no_regions[image]:
        if len(no_region) > 0:
            no_region_low, no_region_hi = no_region
            Ha_backgroundonly[no_region_low:no_region_hi] = np.nan
            NII_backgroundonly[no_region_low:no_region_hi] = np.nan

    T_xpix.append(x_pix)
    T_Ha.append(Ha_backgroundonly)
    T_NII.append(NII_backgroundonly)

T_xpix_transpose = np.array(T_xpix).transpose()
T_Ha_transpose = np.array(T_Ha).transpose()
T_NII_transpose = np.array(T_NII).transpose()


#MODE + MEAN RUNNING BACKGROUNDS
pix_line_all = []
mode_val_all_Ha = []
mode_val_all_NII = []
mean_val_all_Ha = []
mean_val_all_NII = []

for i, [pix_line, Ha_line, NII_line] in enumerate(zip(T_xpix_transpose,T_Ha_transpose,T_NII_transpose)):

    xpix = pix_line[0]
    modeval_Ha = mode(Ha_line)[0][0]
    modeval_NII = mode(NII_line)[0][0]

    pix_line_all.append(xpix)
    mode_val_all_Ha.append(modeval_Ha)
    mode_val_all_NII.append(modeval_NII)
    mean_val_all_Ha.append(np.nanmean(Ha_line))
    mean_val_all_NII.append(np.nanmean(NII_line))


#Backgrounds into Dictionary
Met_N2['Backgrounds']['mean']['Ha'] = mean_val_all_Ha
Met_N2['Backgrounds']['mean']['NII'] = mean_val_all_NII
Met_N2['Backgrounds']['mode']['Ha'] = mode_val_all_Ha
Met_N2['Backgrounds']['mode']['NII'] = mode_val_all_NII


#CALCULATING N2, log(O/H) + 12 (linear and 2nd order fit), and adding to Met_N2 Dictionary
#Building background subtracted Halpha and NII
for i, [image,x_pix, Ha, NII] in enumerate(zip(Met_N2['image'],Met_N2['x_pix'],
                                               Met_N2['Halpha6563*'],Met_N2['[NII]6584weak'])):

    NII_touse = NII - mode_val_all_NII
    Ha_touse = Ha - mode_val_all_Ha
    Met_N2['NII_bkgdsub'].append(NII_touse)
    Met_N2['Halpha_bkgdsub'].append(Ha_touse)

    N2 = func_N2(NII_touse, Ha_touse)
    logOH = func_logOH(N2)
    logOH_lin = func_logOH_linear(N2)

    Met_N2['N2'].append(N2)
    Met_N2['logOH_linear'].append(logOH_lin)
    Met_N2['logOH'].append(logOH)
    Met_N2['HighSignalRegions'].append(dictionary_no_regions[image])
    Met_N2['StarRegions'].append(dictionary_star_regions[image])


#building low-signal backgrounds and sigmas in NII, Halpha
for i, [image, x_pix, Ha, NII] in enumerate(zip(Met_N2['image'],Met_N2['x_pix'],
                                                Met_N2['Halpha_bkgdsub'],Met_N2['NII_bkgdsub'])):
    Ha_backgroundonly = Ha.copy()
    NII_backgroundonly = NII.copy()
    for no_region in dictionary_no_regions[image]:
        if len(no_region) > 0:
            no_region_low, no_region_hi = no_region
            Ha_backgroundonly[no_region_low:no_region_hi] = np.nan
            NII_backgroundonly[no_region_low:no_region_hi] = np.nan

    T_xpix.append(x_pix)
    T_Ha.append(Ha_backgroundonly)
    T_NII.append(NII_backgroundonly)

    #Saving background regions
    Met_N2['Backgrounds']['regions_bkgdsub']['Ha'].append(Ha_backgroundonly)
    Met_N2['Backgrounds']['regions_bkgdsub']['NII'].append(NII_backgroundonly)

    #Sigmas for each line into dictionary
    Met_N2['Backgrounds']['sigma_bkgdsub']['Ha'].append(np.nanstd(Ha_backgroundonly))
    Met_N2['Backgrounds']['sigma_bkgdsub']['NII'].append(np.nanstd(NII_backgroundonly))


#running sigmas - background linearly fit + std (of bkgd regions) above background-sub Ha/NII
def find_lin_fit(X, Y, init_vals = [0, 0, 0]):

    def linear(x, x0, a, b):
        return a * (x - x0) + b

    x_nonan = [a for a,y in enumerate(Y) if math.isnan(y) == False]
    y_nonan = [y for y in Y if math.isnan(y) == False]

    best_vals, covar = curve_fit(linear, x_nonan, y_nonan, p0=init_vals)

    y_fit_vals = linear(X, best_vals[0], best_vals[1], best_vals[2])
    return best_vals, y_fit_vals

Background_Fits_All = {}
Background_Fits_All['Ha'] = []
Background_Fits_All['NII'] = []

for i, [image, x_pix, Ha_bkgdsub, NII_bkgdsub, Ha_bkgdsub_bkgds, NII_bkgdsub_bkgds, Ha_sig, NII_sig] in enumerate(zip(Met_N2['image'],Met_N2['x_pix'],
                                                                                                                      Met_N2['Halpha_bkgdsub'],Met_N2['NII_bkgdsub'],
                                                                                                                      Met_N2['Backgrounds']['regions_bkgdsub']['Ha'],
                                                                                                                      Met_N2['Backgrounds']['regions_bkgdsub']['NII'],
                                                                                                                      Met_N2['Backgrounds']['sigma_bkgdsub']['Ha'],
                                                                                                                      Met_N2['Backgrounds']['sigma_bkgdsub']['NII'])):


    best_vals_x_Ha, y_fit_vals_x_Ha = find_lin_fit(x_pix, Ha_bkgdsub_bkgds)
    best_vals_x_NII, y_fit_vals_x_NII = find_lin_fit(x_pix, NII_bkgdsub_bkgds)

    Background_Fits_All['Ha'].append(best_vals_x_Ha)
    Background_Fits_All['NII'].append(best_vals_x_NII)

    Met_N2['Halpha_sigmas'].append((Ha_bkgdsub - y_fit_vals_x_Ha)/Ha_sig)
    Met_N2['NII_sigmas'].append((NII_bkgdsub - y_fit_vals_x_NII)/NII_sig)

Met_N2['Backgrounds']['background_fits_all'] = Background_Fits_All

# #Running sigmas
# for Ha_bkgdsub, NII_bkgdsub, Ha_sig, NII_sig in zip(Met_N2['Halpha_bkgdsub'],Met_N2['NII_bkgdsub'],
#                                                     Met_N2['Backgrounds']['sigma_bkgdsub']['Ha'],
#                                                     Met_N2['Backgrounds']['sigma_bkgdsub']['NII']):
#     Met_N2['Halpha_sigmas'].append(Ha_bkgdsub/Ha_sig)
#     Met_N2['NII_sigmas'].append(NII_bkgdsub/NII_sig)


#RA DEC + NWSE Dispersion + RA arrays, DEC arrays for all images

#Reading in RADEC files
dir_May31 = '/Users/lockepatton/Desktop/Research/Levesque/NGC6946/datafiles/'
SNcoord = Table.read(dir_May31+'starlist_0530_notes_NGC6946_SNfiles_coord_match.txt',format='ascii')

coord,ra,dec = findRADEC(SNcoord['RAh'],SNcoord['RAm'],SNcoord['RAs'],\
                         SNcoord['DECd'],SNcoord['DECm'],SNcoord['DECs'])


# BUILDING RA DEC in degrees, Dispersion Direction (NS), RA DEC Arrays for All Images
for i,[i1,i2] in enumerate(zip(Met_N2['image'],SNcoord['fitsimage'][::2])):
    Met_N2['RAh'].append(SNcoord['RAh'][::2][i])
    Met_N2['RAm'].append(SNcoord['RAm'][::2][i])
    Met_N2['RAs'].append(SNcoord['RAs'][::2][i])
    Met_N2['DECd'].append(SNcoord['DECd'][::2][i])
    Met_N2['DECm'].append(SNcoord['DECm'][::2][i])
    Met_N2['DECs'].append(SNcoord['DECs'][::2][i])

    RA = ra[::2][i]
    DEC = dec[::2][i]
    Met_N2['RA'].append(RA)
    Met_N2['DEC'].append(DEC)
    Met_N2['NWSE'].append('NS')

    length = 721

    #RA
    RA_array = np.linspace(RA, RA, length)
    #doesn't change

    #DEC
    DEC_low = DEC + 0.05
    DEC_high = DEC - 0.05
    DEC_array = np.linspace(DEC_low, DEC_high, length)
    #changes by -slit_length_u.value/2 DECminutes each way  =  0.05 deg

    Met_N2['RA_array'].append(RA_array)
    Met_N2['DEC_array'].append(DEC_array)



# PLOTS

#Plotting Selected (by hand) High Signal Regions
fig, ax = plt.subplots(len(Met_N2['image']),1)
fig.set_size_inches(10*2,15*1.5)
for i, [image, x_pix, Ha, NII] in enumerate(zip(Met_N2['image'],Met_N2['x_pix'],Met_N2['Halpha6563*'],Met_N2['[NII]6584weak'])):
    ax[i].plot(x_pix, Ha, alpha=.75, c='blue')
    ax[i].plot(x_pix, NII, alpha=.75, c='green')

    ax[i].set_xlabel('Pixel')
    ax[i].set_ylabel(image[8:])
    ax[i].set_title(str(i))
    ax[i].set_ylim(.25e-14, .4e-14)
    ax[i].set_xlim(0,721)

    for no_region in dictionary_no_regions[image]:
        if len(no_region) > 0:
            ax[i].axvspan(no_region[0], no_region[1], alpha=.5, color='k')
fig.savefig('./plots/NGC6946_May31_HighSignalRegions.png')



# Plotting Histograms of all pixel locations low-signal backgrounds + their mode/mean histograms
fig, ax = plt.subplots(1,1)
fig.set_size_inches(10,10)

n_bins = 40

def colormap_f_calc(vals):
    minsub = vals-min(vals)
    return (minsub/max(minsub))
colormap = cm.viridis(colormap_f_calc(x_pix))

for i, [pix_line, Ha_line, NII_line, color] in enumerate(zip(T_xpix_transpose,T_Ha_transpose,T_NII_transpose,colormap)):

    xpix = pix_line[0]
    modeval_Ha = mode(Ha_line)[0][0]
    modeval_NII = mode(NII_line)[0][0]

    Ha_line_nonan = [x for x in Ha_line if str(x) != 'nan']
#     Ha_line_nonan = np.nan_to_num(Ha_line)
    ax.hist(Ha_line_nonan, n_bins, normed=1, histtype='step', stacked=False, fill=False, alpha=.2, color=color)

ax.hist(mode_val_all_Ha, 200, normed=1, histtype='step', stacked=False, fill=False, color='r', label='mode');
ax.hist(mean_val_all_Ha, 200, normed=1, histtype='step', stacked=False, fill=False, color='b', label='mean');

ax.legend();
fig.savefig('./plots/NGC6946_May31_LowSignal_Background_Histograms.png')



#Plotting Running Mode, Mean Backgrounds for NII and Halpha
fig, ax  = plt.subplots()
ax.plot(pix_line_all, mode_val_all_Ha,label='mode - Ha', c='red')
ax.plot(pix_line_all, mean_val_all_Ha,label='mean - Ha', c='darkred')
ax.plot(pix_line_all, mode_val_all_NII,label='mode - NII',c='blue')
ax.plot(pix_line_all, mean_val_all_NII,label='mean - NII',c='darkblue')
ax.legend()
ax.set_ylabel('Running Background')
ax.set_xlabel('Pixel Location on Profile');
fig.savefig('./plots/NGC6946_May31_Running_Backgrounds_Mode+Mean.png')



#plotting final NII, Hapha before + after Background Subtraction
fig, ax = plt.subplots(1,1)
fig.set_size_inches(10*2,4*4)
ax.set_xlabel('Pixel / Galaxy locations')
ax.set_ylabel('Halpha + NII Profiles +offset')

y_baseline = np.linspace(1,1,len(x_pix))

for i, [image,x_pix, Ha, NII] in enumerate(zip(Met_N2['image'],Met_N2['x_pix'],
                                               Met_N2['Halpha6563*'],Met_N2['[NII]6584weak'])):
    baseshift = 1.3e-15

    to_subtract_values = [min(Ha), min(NII)]
    Ha_toplot = Ha - to_subtract_values[0] + y_baseline*i*baseshift
    NII_toplot = NII - to_subtract_values[1] + y_baseline*i*baseshift

    Ha_toplot_running_background = Ha - mode_val_all_Ha + y_baseline*i*baseshift
    NII_toplot_running_background = NII - mode_val_all_NII + y_baseline*i*baseshift

    ax.plot(x_pix, Ha_toplot, alpha=.75, c='blue')
    ax.plot(x_pix, NII_toplot, alpha=.75, c='green')

    ax.plot(x_pix, Ha_toplot_running_background, alpha=.75, c='darkblue')
    ax.plot(x_pix, NII_toplot_running_background, alpha=.75, c='darkgreen')

    ax.text(x_pix[-1]+10, i*baseshift, image)

ax.grid()
fig.savefig('./plots/NGC6946_May31_NII_Ha_After_Mode_Running_Background_Sub.png')



#NII / Halpha Sigma Plots with Background regions + full signals
length_all = len(Met_N2['Halpha_bkgdsub'])
length = 721

fig, ax = plt.subplots(length_all,1)
fig.set_size_inches(20,length_all*3)

i = 0
for Ha_bkgdsub, NII_bkgdsub, Ha_sig, NII_sig, bkgreg_Ha, bkgreg_NII, x_pix, im in zip(Met_N2['Halpha_bkgdsub'],Met_N2['NII_bkgdsub'],
                                                                                      Met_N2['Backgrounds']['sigma_bkgdsub']['Ha'],
                                                                                      Met_N2['Backgrounds']['sigma_bkgdsub']['NII'],
                                                                                      Met_N2['Backgrounds']['regions_bkgdsub']['Ha'],
                                                                                      Met_N2['Backgrounds']['regions_bkgdsub']['NII'],
                                                                                      Met_N2['x_pix'], Met_N2['image']):
    ax[i].plot(x_pix,Ha_bkgdsub, 'b')
    ax[i].plot(bkgreg_Ha, 'r', alpha=.5)
    ax[i].plot(x_pix,NII_bkgdsub, 'g')
    ax[i].plot(bkgreg_NII, 'r', alpha=.5)

    for sig_scale_num in range(5):
        ax[i].plot(x_pix,[Ha_sig*sig_scale_num]*721, c='darkblue', alpha=0.4)
        ax[i].plot(x_pix,[NII_sig*sig_scale_num]*721,  c='darkgreen', alpha=0.4)

    ax[i].set_xlim(0,721)
    ax[i].set_ylim(-Ha_sig, 7*Ha_sig)
    ax[i].set_title(str(i) +' | '+im)
    ax[i].grid()

    i+=1

fig.savefig('./plots/NGC6946_May31_Halpha_NII_Sigma_Plots.png')


#SAVING DICTIONARY
np.save('Dictionary_y_NGC6946_May31.npy', Met_N2)

__author__ = 'Locke Patton'
# Run 'ulimit -n 512' before executing

# CHANNELS - EMISSION LINES

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

# Directories
dir_5apally = '/Users/lockepatton/Desktop/Research/Levesque/NG6946_May24/data/attempt2/5apally/'
dir_trace_centers = '/Users/lockepatton/Desktop/Research/Levesque/NG6946_May24/jupyter/apall_y/'
#Reading in y apall Trace Extract Centers from TraceCenters_All.xlsx
Excel_Centers = pd.read_excel(dir_trace_centers+'/TraceCenters_All.xlsx')


#function that creates a SlitSpectra Object in y direction
def CreateSS_y(image, trace_center=0, plotSpec=False, verbose=False):
    channel = image[-1]
    if channel == 'r':
        dir_ = dir_5apally + 'red/'
    elif channel == 'b':
        dir_ = dir_5apally + 'blue/'
    full_region = [1,2098,1,1078]
    center_width = 10

    trim_region_r=[15,1920,99,819]
    full_lowest_r=83.5
    full_highest_r=834.5

    trim_region_b=[210,1990,210,888]
    full_lowest_b=195
    full_highest_b=903.5

    if channel == 'r':
#         trace_center = 65.29
        trim_region = trim_region_r
        full_lowest = full_lowest_r
        full_highest = full_highest_r
#         w1 =  5428.431 #starting wavelength
#         w2 =  9802.673 #ending wavelength
#         dw =  2.29619 #wavelength intervel per pixel
#         nw =  1906 #number of output pixels
        w1 =  5461.809 #starting wavelength  #May 24th Night
        w2 =  9837.645 #ending wavelength
        dw =  2.297027 #wavelength intervel per pixel
        nw =  1906 #number of output pixels

    if channel == 'b':
#         trace_center = 1458.02
        trim_region = trim_region_b
        full_lowest = full_lowest_b
        full_highest = full_highest_b
#         w1 =  2926.267 #starting wavelength
#         w2 =  6165.38 #ending wavelength
#         dw =  1.819726 #wavelength intervel per pixel
#         nw =  1781 #number of output pixels
        w1 =  2906.903 #starting wavelength  #May 24th Night
        w2 =  6157.6 #ending wavelength
        dw =  1.826235 #wavelength intervel per pixel
        nw =  1781 #number of output pixels

#     trim_base_vs_transfer_fit = [0.94335225214255913, 0.81555279062263253]
#     c1,c2 = trim_base_vs_transfer_fit

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
        T.plotImageWavelength(save_fig=True, verbose=False, blackout=False,vmin=0, vmax=1.7e-15)
        T.plotImageWavelength(save_fig=True, verbose=False, blackout=False,vmin=0, vmax=1.7e-15, x=[300,600])

    T.buildSpectra(verbose=verbose)

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


for image, center in zip(Excel_Centers['File'],Excel_Centers['Trace Center']):

    image = image[:-5]
    T = CreateSS_y(image,trace_center=center,plotSpec=False)

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
    'fdFT4.1_1917A_EW.0019r' : [[110,160], [400,620], [680,700]],
    'fdFT4.1_1917A_NS.0018r' : [[0,150], [290,500]],
    'fdFT4.1_1917A_NS.0040r' : [[0,140], [340,400], [530,570], [660,700]],
    'fdFT4.1_1939C_EW.0025r' : [[0,60], [120,130], [330,350], [700,721]],
    'fdFT4.1_1939C_NS.0024r' : [[10,50], [50,60], [250,280], [310,320], [510,560]],
    'fdFT4.1_1939C_NS.0037r' : [[50,60], [260,320], [580,600]],
    'fdFT4.1_1948B_EW.0013r' : [[250,340], [430,560]],
    'fdFT4.1_1948B_NS.0012r' : [[340,560], [610,630]],
    'fdFT4.1_1948B_NS.0035r' : [[340,370], [410,550], [610,630]],
    'fdFT4.1_1968D_EW.0015r' : [[50,80], [140,155], [190,240], [260,280], [300,400], [440,450], [540,560], [580,720]],
    'fdFT4.1_1968D_NS.0014r' : [[210,300], [515,530], [590,640], [680,700]], #*
    'fdFT4.1_1968D_NS.0042r' : [[120,130], [230,260], [280,310], [510,615], [660,680]],
    'fdFT4.1_1969P_EW.0017r' : [[120,130], [265,290], [450,470]],
    'fdFT4.1_1969P_NS.0016r' : [[100,140], [160,200], [240,260], [310,350], [435,445], [455,465], [530,560], [660,680]],
    'fdFT4.1_1969P_NS.0041r' : [[0,40], [85,105], [220,270], [280,310], [500,560]],
    'fdFT4.1_1980K_EW.0006r' : [[160,180]],
    'fdFT4.1_1980K_NS.0005r' : [[0,10]],
    'fdFT4.1_1980K_NS.0032r' : [[30,65], [120,130]],
    'fdFT4.1_1980K_NS.0034r' : [[260,290], [430,450], [700,721]],
    'fdFT4.1_2002hh_EW.0023r' : [[230,255], [320,340], [400,570], [610,700]],
    'fdFT4.1_2002hh_NS.0022r' : [[420,450], [700,721]],
    'fdFT4.1_2002hh_NS.0038r' : [[20,80], [115,135], [170,200], [420,470]],
    'fdFT4.1_2004et_EW.0010r' : [[100,140], [240,250]],
    'fdFT4.1_2004et_NS.0009r' : [[0,100]],
    'fdFT4.1_2004et_NS.0033r' : [[0,60], [110,135], [330,355]],
    'fdFT4.1_2008S_EW.0021r' : [[20,50], [230,250]],
    'fdFT4.1_2008S_NS.0020r' : [[200,280], [620,640]],
    'fdFT4.1_2008S_NS.0039r' : [[70,90], [210,240], [700,721]],
    'fdFT4.1_L5_NS.0047r' : [[120,130]],
    'fdFT4.1_L6_NS.0048r' : [[60,80], [190,270], [520,560]],
    'fdFT4.1_L7_NS.0049r' : [[50,60], [120,130], [170,190], [290,330], [430,450]],
    'fdFT4.1_L8_NS.0050r' : [[0,40]],
    'fdFT4.1_LBV_NS.0026r' : [[220,340], [490,510]],
    'fdFT4.1_LBV_NS.0036r' : [[100,140], [220,310], [420,450], [530,560], [660,700]],
    'fdFT4.1_MF16_EW.0031r' : [[0,40], [430,470], [590,670]],
    'fdFT4.1_MF16_NS.0030r' : [[490,721]],
    'fdFT4.1_U5_NS.0043r' : [[0,80], [580,721]],
    'fdFT4.1_U6_NS.0044r' : [[250,280], [335,370], [580,721]],
    'fdFT4.1_U7_NS.0045r' : [[120,130], [340,390], [435,450], [610,721]],
    'fdFT4.1_U8_NS.0046r' : [[190,210], [640,680]],
    'fdFT4.1_center_EW.0029r' : [[50,60], [80,210], [260,300], [380,550], [590,640]],
    'fdFT4.1_center_NS.0027r' : [[50,60], [190,450], [530,580], [630,721]],
    'fdFT4.1_center_NS.0028r' : [[50,60], [260,290], [310,320], [350,490], [530,560], [650,721]],
}

#Stars + Bad Pixels

dictionary_star_regions = {
    'fdFT4.1_1917A_EW.0019r' : [[135,150]], #might be a star - check - NII ~ Halpha
    'fdFT4.1_1917A_NS.0018r' : [],
    'fdFT4.1_1917A_NS.0040r' : [],
    'fdFT4.1_1939C_EW.0025r' : [],
    'fdFT4.1_1939C_NS.0024r' : [[12,30]], #might be a star - check - NII ~ Halpha
    'fdFT4.1_1939C_NS.0037r' : [],
    'fdFT4.1_1948B_EW.0013r' : [],
    'fdFT4.1_1948B_NS.0012r' : [],
    'fdFT4.1_1948B_NS.0035r' : [],
    'fdFT4.1_1968D_EW.0015r' : [[447,452]], #bad pixel
    'fdFT4.1_1968D_NS.0014r' : [[515,529]], #might be a star - check - NII ~ Halpha
    'fdFT4.1_1968D_NS.0042r' : [[660,680]], #might be a star - check - NII ~ Halpha
    'fdFT4.1_1969P_EW.0017r' : [],
    'fdFT4.1_1969P_NS.0016r' : [],
    'fdFT4.1_1969P_NS.0041r' : [[89,100]], #star
    'fdFT4.1_1980K_EW.0006r' : [],
    'fdFT4.1_1980K_NS.0005r' : [],
    'fdFT4.1_1980K_NS.0032r' : [],
    'fdFT4.1_1980K_NS.0034r' : [[265,287]],
    'fdFT4.1_2002hh_EW.0023r' : [],
    'fdFT4.1_2002hh_NS.0022r' : [[420,440]], #star??? #HALP Don't know whether these are stars
    'fdFT4.1_2002hh_NS.0038r' : [[433,453]], #[[430,455]], Probably not a star? NII < Halpha, but less than in other emission regions
    'fdFT4.1_2004et_EW.0010r' : [],
    'fdFT4.1_2004et_NS.0009r' : [],
    'fdFT4.1_2004et_NS.0033r' : [],
    'fdFT4.1_2008S_EW.0021r' : [[25,45]], #def a star
    'fdFT4.1_2008S_NS.0020r' : [],
    'fdFT4.1_2008S_NS.0039r' : [[700,721]], #def a star, but in a region with multiple bumps - leave 2nd bump?
    'fdFT4.1_L5_NS.0047r' : [],
    'fdFT4.1_L6_NS.0048r' : [[64,77]], #star
    'fdFT4.1_L7_NS.0049r' : [],
    'fdFT4.1_L8_NS.0050r' : [],
    'fdFT4.1_LBV_NS.0026r' : [],
    'fdFT4.1_LBV_NS.0036r' : [[428,450]], #def a star
    'fdFT4.1_MF16_EW.0031r' : [[450,465]], #probably a star? but NII changes across profile
    'fdFT4.1_MF16_NS.0030r' : [],
    'fdFT4.1_U5_NS.0043r' : [[50,75],[700,715]], #1st def yes #2nd prob not
    'fdFT4.1_U6_NS.0044r' : [[250,272]], #def a star
    'fdFT4.1_U7_NS.0045r' : [],
    'fdFT4.1_U8_NS.0046r' : [],
    'fdFT4.1_center_EW.0029r' : [],
    'fdFT4.1_center_NS.0027r' : [],
    'fdFT4.1_center_NS.0028r' : [[537,552]], #Star within HII region? (check with Emily)
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

# #Running sigmas (from 0)
# for Ha_bkgdsub, NII_bkgdsub, Ha_sig, NII_sig in zip(Met_N2['Halpha_bkgdsub'],Met_N2['NII_bkgdsub'],
#                                                     Met_N2['Backgrounds']['sigma_bkgdsub']['Ha'],
#                                                     Met_N2['Backgrounds']['sigma_bkgdsub']['NII']):
#     Met_N2['Halpha_sigmas'].append(Ha_bkgdsub/Ha_sig)
#     Met_N2['NII_sigmas'].append(NII_bkgdsub/NII_sig)


#RA DEC + NWSE Dispersion + RA arrays, DEC arrays for all images

#Reading in RADEC files
dir_May24 = '/Users/lockepatton/Desktop/Research/Levesque/NG6946_May24/data/attempt2/3wcalibration/'
RADEC_May24 = pd.read_csv(dir_May24 + 'RADECNGC6946_May24.txt', delimiter='\t')

# BUILDING RA DEC in degrees
for index in RADEC_May24.index:
    ra_,dec_ = RADEC_May24['RA'][index],RADEC_May24['DEC'][index]
    RAh, RAm, RAs = ra_.strip().split(':')
    DECd, DECm, DECs = dec_.strip().split(':')
    Met_N2['RAh'].append(int(RAh))
    Met_N2['RAm'].append(int(RAm))
    Met_N2['RAs'].append(float(RAs))
    Met_N2['DECd'].append(int(DECd))
    Met_N2['DECm'].append(int(DECm))
    Met_N2['DECs'].append(float(DECs))

RADECMay24_coord, RADECMay24_ra, RADECMay24_dec = findRADEC(Met_N2['RAh'],Met_N2['RAm'],Met_N2['RAs'],
                                                            Met_N2['DECd'],Met_N2['DECm'],Met_N2['DECs'])

Met_N2['RA'] = RADECMay24_ra
Met_N2['DEC'] = RADECMay24_dec

#Determining Dispersion Direction (NS or EW)
def determineEWNS(angle):
    if angle == 0.0:
        return 'EW'
    if angle == 90.0:
        return 'NS'
    else:
        print angle,'- EWNS of angle was not determined.'

Met_N2['NWSE'] = map(determineEWNS,RADEC_May24['OBJANGLE'])

#Building RA DEC Arrays for All Images
for NWSE,RA,DEC in zip(Met_N2['NWSE'],Met_N2['RA'],Met_N2['DEC']):
    length = 721
    if NWSE == 'NS':
        #RA
        RA_array = np.linspace(RA, RA, length)
        #doesn't change

        #DEC
        DEC_low = DEC + 0.05
        DEC_high = DEC - 0.05
        DEC_array = np.linspace(DEC_low, DEC_high, length)
        #changes by -slit_length_u.value/2 DECminutes each way  =  0.05 deg
    elif NWSE == 'EW':
        #RA
        RA_low = RA - 0.05
        RA_high = RA + 0.05
        RA_array = np.linspace(RA_low, RA_high, length)
        #changes by +24./2. RAseconds each way  = 0.05 deg

        #DEC
        DEC_array = np.linspace(DEC, DEC, length)
        #doesn't change
    else:
        print 'error : no calculated RA/DEC'
        RA_array=[]
        DEC_array=[]

    Met_N2['RA_array'].append(RA_array)
    Met_N2['DEC_array'].append(DEC_array)




#PLOTS

#Plotting Selected (by hand) High Signal Regions
fig, ax = plt.subplots(len(Met_N2['image']),1)
fig.set_size_inches(10*2,8*10*1.5)
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
fig.savefig('./plots/NGC6946_May24_HighSignalRegions.png')



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
fig.savefig('./plots/NGC6946_May24_LowSignal_Background_Histograms.png')



#Plotting Running Mode, Mean Backgrounds for NII and Halpha
fig, ax  = plt.subplots()
ax.plot(pix_line_all, mode_val_all_Ha,label='mode - Ha', c='red')
ax.plot(pix_line_all, mean_val_all_Ha,label='mean - Ha', c='darkred')
ax.plot(pix_line_all, mode_val_all_NII,label='mode - NII',c='blue')
ax.plot(pix_line_all, mean_val_all_NII,label='mean - NII',c='darkblue')
ax.legend()
ax.set_ylabel('Running Background')
ax.set_xlabel('Pixel Location on Profile');
fig.savefig('./plots/NGC6946_May24_Running_Backgrounds_Mode+Mean.png')



#plotting final NII, Hapha before + after Background Subtraction
fig, ax = plt.subplots(1,1)
fig.set_size_inches(10*2,8*4)
ax.set_xlabel('Pixel / Galaxy Locations')
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
fig.savefig('./plots/NGC6946_May24_NII_Ha_After_Mode_Running_Background_Sub.png')



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

fig.savefig('./plots/NGC6946_May24_Halpha_NII_Sigma_Plots.png')


#SAVING DICTIONARY
np.save('Dictionary_y_NGC6946_May24.npy', Met_N2)

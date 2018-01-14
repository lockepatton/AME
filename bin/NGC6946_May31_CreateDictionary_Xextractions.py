from SlitSpectra4 import SlitSpectra as SS
import warnings
import os
import numpy as np
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
from scipy.stats import mode


# Verbose
verbose = False
verbose_absolute = True

# Directories
# line profile inputs
# location of slot.log files
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


#reading in apall line fit data splot.logs and separating into new files based upon aperture.
dictionary_of_splotlogs = {}
dictionary_of_splotlogs['image'] = []
dictionary_of_splotlogs['table'] = []

for cdFTImage_r_ in cdFTImages_r:

    splot_log = dir_lineprofile_ + cdFTImage_r_ + '.splot.log'

    if os.path.isfile(splot_log):
        with open(splot_log, 'r') as file:

            data_splot_log = file.readlines()

            starts = []
            for i, base in enumerate(data_splot_log):
                if len([letterindex for letterindex, letter in enumerate(base) if letter == '[']) != 0:
                    starts.append(i)

            #making splotlog directory
            directory = dir_lineprofile_ + 'splotlog/'
            if not os.path.exists(directory):
                os.makedirs(directory)

            #making indecies simply for index of "data", wich includes the last index as len(data_splot_log)
            indecies = starts[:]
            indecies.append(len(data_splot_log))

            #saving all splot.log separated by aperture
            for i_,start_ in enumerate(starts):

                sub_sub_image_name = str(data_splot_log[start_][:-3].replace(' ','_'))
                sub_image_name = str(cdFTImage_r_) + '_' + sub_sub_image_name

                to_save = dir_lineprofile_ + 'splotlog/' + sub_image_name

                if verbose:
                    print to_save

                #writing individual files with center, cont, flux, eqw, core, gftwhm, lfwhm values
                with open(to_save,'w') as file:
                    data = data_splot_log[start_:indecies[i_+1]]
                    data[0] = '   center      cont      flux       eqw      core     gfwhm     lfwhm \n'
                    file.writelines(data)

                #reading in same files as astropy.table
                T = Table.read(to_save,format='ascii')

                #adding said table and image name to dictionary for later use :D
                dictionary_of_splotlogs['image'].append(sub_image_name)
                dictionary_of_splotlogs['table'].append(T)





#BUILDING Dictionary separated by Images

dictionary_of_ap_N2_metallicity_byimage = {}

#building byimage dictionary
for im in cdFTImages_r:
    dictionary_of_ap_N2_metallicity_byimage[im]={}
    for name_ in ['savename','aperture','Ha','NII','N2','logOH','logOH_linear']:
        dictionary_of_ap_N2_metallicity_byimage[im][name_]=[]

Halpha = 6563
NII = 6584
deltaHaNII = Halpha - NII

func_N2 = lambda NII,Halpha : np.log10(NII/Halpha)
func_logOH = lambda N2 : 9.37 + 2.03*N2 + 1.26*N2**2 + 0.32*N2**3
func_logOH_linear = lambda N2 : 8.90 + 0.57* N2

if verbose:
    print 'true Ha - NII', deltaHaNII

#TODO build this better
#i.e. instead of using tolerance alone, find all delta_Ha_NII values and look for closest within range
Angstrom_err_tolerance = 2

print 'aperture, image, savename, Ha-NII'
print 'NII[center,flux], Halpha[center,flux], N2, logOH, logOH_linear','\n'

for save_name,test_kk in zip(dictionary_of_splotlogs['image'],dictionary_of_splotlogs['table']):

    #TODO build this better
    #sorting via flux. assuming 2nd NII and Halpha have highest flux - Replace
    test_kk.sort('flux')

    #defining im (image name) from inside save_name
    mark_locations = [letterindex for letterindex, letter in enumerate(save_name) if letter == '_']
    im = save_name[:mark_locations[2]]

    if (len(test_kk)>=2):
        Halpha_candidate = test_kk[-1]
        NII_candidate =  test_kk[-2]
        deltaHaNII_candidate = Halpha_candidate['center']-NII_candidate['center']

        if verbose:
            print 'candidate:', deltaHaNII_candidate

        if ((deltaHaNII_candidate-deltaHaNII)**2 < Angstrom_err_tolerance**2):

            N2 = func_N2(NII_candidate['flux'],Halpha_candidate['flux'])

            dictionary_of_ap_N2_metallicity_byimage[im]['savename'].append(save_name)
            dictionary_of_ap_N2_metallicity_byimage[im]['aperture'].append(save_name[-13-4:-13])
            dictionary_of_ap_N2_metallicity_byimage[im]['Ha'].append(Halpha_candidate)
            dictionary_of_ap_N2_metallicity_byimage[im]['NII'].append(NII_candidate)
            dictionary_of_ap_N2_metallicity_byimage[im]['N2'].append(N2)
            dictionary_of_ap_N2_metallicity_byimage[im]['logOH'].append(func_logOH(N2))
            dictionary_of_ap_N2_metallicity_byimage[im]['logOH_linear'].append(func_logOH_linear(N2))

            if verbose_absolute:
                print save_name[-13-4:-13],save_name[:mark_locations[2]], save_name, deltaHaNII_candidate
                print '[',NII_candidate['center'],NII_candidate['flux'],'] [', \
                      Halpha_candidate['center'],Halpha_candidate['flux'],']', \
                      N2,func_logOH(N2),func_logOH_linear(N2),'\n'

if verbose_absolute:
    for im in cdFTImages_r:
        print im, len(dictionary_of_ap_N2_metallicity_byimage[im]['aperture']),\
              dictionary_of_ap_N2_metallicity_byimage[im]['aperture']


#SAVING DICTIONARY
np.save('Dictionary_x_NGC6946_May31.npy', dictionary_of_ap_N2_metallicity_byimage)

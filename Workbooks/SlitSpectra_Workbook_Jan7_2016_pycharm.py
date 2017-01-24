from SlitSpectra2 import SlitSpectra2 as SS
import warnings
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
warnings.filterwarnings('ignore')
#%pylab
#%matplotlib inline

plt.plot([1,2,3,4],[1,4,9,16], 'ro')

dTImages_r = ['dT4.1_02hh_NS.0052r',
            'dT4.1_04et_NS.0047r',
            'dT4.1_08S_NS.0051r',
            'dT4.1_17A_NS.0050r',
            'dT4.1_39C_NS.0053r',
            # 'dT4.1_48B_NS.0046r',
            'dT4.1_68D_NS.0045r',
            'dT4.1_69P_NS.0049r',
            'dT4.1_80K_NS.0048r']

dFTImages_r = ['dFT4.1_02hh_NS.0052r',
            'dFT4.1_04et_NS.0047r',
            'dFT4.1_08S_NS.0051r',
            'dFT4.1_17A_NS.0050r',
            'dFT4.1_39C_NS.0053r',
            # 'dFT4.1_48B_NS.0046r',
            'dFT4.1_68D_NS.0045r',
            'dFT4.1_69P_NS.0049r',
            'dFT4.1_80K_NS.0048r']

dTImages_b = ['dT4.1_02hh_NS.0052b',
            'dT4.1_04et_NS.0047b',
            'dT4.1_08S_NS.0051b',
            'dT4.1_17A_NS.0050b',
            'dT4.1_39C_NS.0053b',
            # 'dT4.1_48B_NS.0046b',
            'dT4.1_68D_NS.0045b',
            'dT4.1_69P_NS.0049b',
            'dT4.1_80K_NS.0048b']

dFTImages_b = ['dFT4.1_02hh_NS.0052b',
            'dFT4.1_04et_NS.0047b',
            'dFT4.1_08S_NS.0051b',
            'dFT4.1_17A_NS.0050b',
            'dFT4.1_39C_NS.0053b',
            # 'dFT4.1_48B_NS.0046b',
            'dFT4.1_68D_NS.0045b',
            'dFT4.1_69P_NS.0049b',
            'dFT4.1_80K_NS.0048b']

for dT_image_r, dFT_image_r, dT_image_b, dFT_image_b in zip(dTImages_r, dFTImages_r, dTImages_b, dFTImages_b):
    #     image = dT_image_r
    image = dT_image_b
    #     image = dFT_image_r
    #     image = dFT_image_b

    channel = image[-1]

    #     image_to_write = dT_image_b

    dir_ = '/Users/amielpatton-hall/Desktop/Research/Levesque/NGC6946/data2.0/apall4/'
    full_region = [1, 2098, 1, 1078]
    center_width = 10
    trace_center = 0  # initially set before know value.

    trim_region_r = [15, 1920, 99, 819]
    full_lowest_r = 83.5
    full_highest_r = 834.5

    trim_region_b = [210, 1990, 210, 888]
    full_lowest_b = 195
    full_highest_b = 903.5

    if channel == 'r':
        trim_region = trim_region_r
        full_lowest = full_lowest_r
        full_highest = full_highest_r
        w1 = 5428.431  # starting wavelength
        w2 = 9802.673  # ending wavelength
        dw = 2.29619  # wavelength intervel per pixel
        nw = 1906  # number of output pixels

    if channel == 'b':
        trim_region = trim_region_b
        full_lowest = full_lowest_b
        full_highest = full_highest_b
        w1 = 2926.267  # starting wavelength
        w2 = 6165.38  # ending wavelength
        dw = 1.819726  # wavelength intervel per pixel
        nw = 1781  # number of output pixels

    if image[-17:] == '4.1_02hh_NS.0052r':
        trace_center = 489
    if image[-17:] == '4.1_04et_NS.0047r':
        trace_center = 159.12  # 157.775
    if image[-16:] == '4.1_08S_NS.0051r':
        trace_center = 712.207
    if image[-16:] == '4.1_17A_NS.0050r':
        trace_center = 313.38
    if image[-16:] == '4.1_39C_NS.0053r':
        trace_center = 652.623
    if image[-16:] == '4.1_48B_NS.0046r':
        trace_center = 501.5273  # 484
    if image[-16:] == '4.1_68D_NS.0045r':
        trace_center = 288.14  # 287.669
    if image[-16:] == '4.1_69P_NS.0049r':
        trace_center = 427.26  # 426.928
    if image[-16:] == '4.1_80K_NS.0048r':
        trace_center = 487.76  # 488.154

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
           w1=w1, w2=w2, dw=dw)

    T.findcenterwidthBackground(verbose=True)
    T.allcalc(verbose=True)

    #     T.plotImage(save_fig=True,verbose=False)
    #     if channel == 'b':
    #         T.plotImage(x=[1000,1200],save_fig=True,verbose=False)  #b Signal in b
    #     if channel == 'r':
    #         T.plotImage(x=[330,700],save_fig=True,verbose=False)  #r Signal in r - Halpha Group

    #     T.buildotherdatabase(image_to_write=image_to_write,verbose=True)

    values_in_r = [125.11111, 255.66667, 344.55556]
    values_in_b = [234.21605, 357.67284, 440.62037]

    #     T.buildothercoordinatefile(transfer_image_to_write=image_to_write,
    #                              base_channel_points=values_in_r,
    #                              transfer_channel_points=values_in_b,
    #                              transfer_full_lowest=full_lowest_b,
    #                              transfer_full_highest=full_highest_b,
    #                              transfer_trim_region=trim_region_b,
    #                              verbose=True)

    #     T.findcenterwidthBackground(verbose=True)
    #     T.allcalc(verbose=True)

    #     T.builddatabase(verbose=True)

    T.plotImage(save_fig=True, verbose=False)
    if channel == 'b':
        T.plotImage(x=[1000, 1200], save_fig=True, verbose=False)  # b Signal in b
    if channel == 'r':
        T.plotImage(x=[330, 700], save_fig=True, verbose=False)  # r Signal in r - Halpha Group


    #     T.buildbasedatabase(verbose=True, database_source=None)    #careful
    #     T.builddatabase(verbose=True)      #careful (not completely fool-proof)

    #     T.plotSpectra(x=[6000,7000],save_fig=False,verbose=False)
    #     T.plotSpectra(verbose=False)

    #     T.finalizeFiles(t='2016-12-28_04h',verbose=False)


    #     print ''
    #     !cat $dir_/database/\ap$image |wc
    #     !cat $dir_/database/\ap$image
    #     print '' ; print ''
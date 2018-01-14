===========
AME  ApallMultislitExtension
===========

This code is designed for use while extracting calibrated spectra with multiple apertures inside IRAF/PyRAF's apall package inside twodspec.apextract. It ultimately takes in a coord file and outputs database files for use inside apall. Also can input wavelengths, and output database files.


There are two main sides to this code:
1. mass extraction help in both x and y dispersion directions when working with apall in IRAF/PyRAF.
2. working with extracted spectra - plotting, reading in mass spectra into python and building data structures, particularly when working with line profiles, but also for spectra.


Example setup
=========
Typical set-up often looks like this::

        #!/usr/bin/env python
        from AME import SlitSpectra as SS

        image = "filename"      #image file name
        dir_ = 'path/to/file/'  #path to image name
        full_region = [1,2098,1,1078]n  #image region dimensions of image before trimming
        center_width = 10       #width of center extraction (code automatically creates center extraction)

        trim_region = [15,1920,99,819]  #trimming region dimensions of image
        full_lowest = 83.5      #lowest pixel in y that data is found
        full_highest = 834.5    #highest pixel in y that data is found

        w1 =  2906.903          #starting wavelength  #May 24th Night
        w2 =  6157.6            #ending wavelength
        dw =  1.826235          #wavelength interval per pixel
        nw =  1781              #number of output pixels

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

Once the SlitSpectra (SS) object is built, there are *many* ways to go from there.



Example Complete Extraction on Images with PyRAF's apall
=========
###### Step 1 - Determine Apertures + Setup Basic Database:

###### In pyraf:

- Create image + \_coords.txt file with upper ap, lower ap, followed by 4 background coords, for each aperture in order from 1+. This can be done inside ds9, using the 'c' key. Order of background coordinates does not matter, but do upper, then lower coordinates in that order.

- Run apall, add one test aperture and save to database, creating a mock base file inside database. At this point, you can change any fundamental components of your extraction. This will generated a file inside your database folder of the form ap+image.

###### In python:
- Define image parameters

###### Step 2 - Find true trace_center:
###### In python: find initial trace_center
- Define an initial trace_center.
    - Run previous cell again until all centers correct.

###### Once satisfied with initial trace_center:
- Add buildbasedatabase() to previous cell to create database file with initial center on 1st aperture

###### In pyraf: use apall to find true trace center + generate good trace function
- Run apall on 1st database file just created.
- Center aperture
- Input this center into python trace_center variable.
- Run trace, save to database, do not extract.

###### Step 3 - In python:
- Run builddatabase()
    - Creates database with all apertures, with correct fitted traces + centers

###### Step 4 - In pyraf:
- Run apall on all apertures
- Check backgrounds and centers
- Extract all apertures, basename = imagename.
- If you made any changes inside apall, make sure to adjust them in your files (since this program won't know).

###### Step 5 - aperture & spectra plots.

- plotSpectra() plotImage() & plotImageWavelength() & plotWavelengthSpectra() can be used to see your data quickly, especially when working with many files at once

###### Step 6 finalizeFiles()
- Finalize files will copy all related files and add any necessary info (including date, time, profile names, etc) to the file names. This way you know exactly when and where and what your files where when you come back to them.
- Make sure to copy your plots into the dated finalized directory, if you want them there too.


Typical setup to run this process looks like this::


Once your profiles + spectra have been extracted
=========


Plotting images
---------

The main command to plot images is plotImageWavelength::        
        T.plotImageWavelength(save_fig=True,verbose=False)
        T.plotImageWavelength(x=[250,700],save_fig=True,verbose=False)

All line profiles that have been defined using T.calculateWavelengths will be added to the plot too, so you can check alignment between where you are extracting (and if you're using them their background regions).


Plotting line profiles
---------

It's easy to plot line profiles::
        T.plotSpectra(save_fig=False, verbose=False)

Typical full usage often looks like this - note that line profiles must have been extracted ::

      wavelengths = [6300, 6548, 6563, 6584, 6717, 6731]     #red wavelength values in Angstroms
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


Returning your spectra / line profile data
---------

Once you've set-up your object and run T.buildSpectra, you can return line profiles / spectra. This process is shown here::

    returnsAllSpec = T.returnAllSpectra()
    [x_pix, SpectraAll] = returnsAllSpec

    Halpha_profile = SpectraAll['Halpha6563*']
    NII_profile = SpectraAll['[NII]6584weak']  #etc


If you'd like to return min-subtracted data, you can also pull that functionality out of T.plotWavelengthSpectra::

    returns = T.plotWavelengthSpectra('[NII]6584weak','Halpha6563*',x=None, minSubtract=False,
                                      verbose=False,returned=True,save_fig=True)

    [x_pix, [y_NII6584_min_sub, y_NII6584_true, wavelength_NII6584], \
            [y_Halpha6563_min_sub, y_Halpha6563_true, wavelength_Halpha6563]] = returns

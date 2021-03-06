from distutils.core import setup

setup(
    name='AME',
    version='0.1.0',
    author='Locke Patton',
    author_email='amielp@uw.edu',
    packages=['AME'],
    scripts=['bin/NGC6946_May24_CreateDictionary.py','bin/NGC6946_May31_CreateDictionary_Xextractions.py', 'bin/NGC6946_May31_CreateDictionary.py'],
    url='https://github.com/lockepatton/AME',
    license='LICENSE.txt',
    description='ApallMultislitExtension is designed for use while extracting calibrated spectra with multiple apertures inside IRAF/PyRAF\'s apall package inside twodspec.apextract, and for further data manipulation of large spectra and line spectra projects.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy",
        "matplotlib",
        "pandas"
    ],
    classifiers=[
        "Programming Language :: Python :: 2 :: Only",
        "License :: MIT License",
        "Intended Audience :: Science/Research",
    ]
)

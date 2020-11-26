# IGMabs

IGMabs: a high resolution X-ray absorption model for highly ionized species (see, for example, [Gatuzz,Garcia and Kallman (2020)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.483L..75G/abstract)). This model is based on ISMabs ([Gatuzz et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...800...29G/abstract))

IGMabs (Intergalactic Medium Absorption) is a high-resolution X-ray photoabsorption model that takes into account the following highly ionized species: CIV-CVI, NV-NVII, OVI-OVIII, NeVIII-NeX, MgX-MgXII, SiXII-SiXIV, SXIV-SXVI, FeXVIII-FeXX, FeXXIV-FeXXVI. The ionic column densities correspond to the free parameters of the model. The model also includes turbulence broadening, which is applied to the optical depth (i.e., all ionic species have the same turbulence velocity)

Particularly we include in our model the following cross sections:

- CIV, CV and CVI from [Hasoglu et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010ApJ...724.1296H/abstract)
- NV, NVI and NVII from [Garcia et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJS..185..477G/abstract)
- OVI, OVII and OVIII from [Garcia et al. (2005)](https://ui.adsabs.harvard.edu/abs/2005ApJS..158...68G/abstract) 
- NeVIII, NeIX and NeX from [Juett et al. (2006)](https://ui.adsabs.harvard.edu/abs/2006ApJ...648.1066J/abstract)
- MgX, MgXI and MgXII from [Hasoglu et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJS..214....8H/abstract)
- SiXII, SiXIII and SiXIV from [Witthoeft et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJS..182..127W/abstract)
- SXIV, SXV and SXVI from [Witthoeft et al. (2009)](https://ui.adsabs.harvard.edu/abs/2009ApJS..182..127W/abstract)
- FeXVIII, FeXIX, FeXX, FeXXIV, FeXXV, FeXXVI from [Bautista et al. (2004)](https://ui.adsabs.harvard.edu/abs/2004A%26A...418.1171B/abstract)
 

OBTAINING IGMabs

The model can be downloaded from the Github repository at https://github.com/efraingatuzz/IGMabs. Please note thaat, because of the large files for the atomic data, you must clone the git repository in your local machine (i.e., if you download a zip file the atomic data will not be included). In order to clone the repository first you must install git and lfs

To install git:

Follow the instructions in https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

To install lfs

1- Navigate through https://git-lfs.github.com/ and click in the download option

2- On your computer, locate and unzip the downloaded file.

3- Open Terminal and change the current working directory into the folder you downloaded and unzipped.

4- To install the file, run the command "./install.sh" (without quotes)

5- Verify that the installation was successful by running the command "git lfs install" (without quotes)

To clone the ISMabs respository:

1- Create a folder in your local machine where the model will be located

2- Run the following command:
git clone https://github.com/efraingatuzz/ISMabs

The contents of the folder include:

- atomic_data/AtomicData.fits  -- atomic database binary fits file. This must reside in the directory atomic_data inside the folder where the model is located.  
- igmabs.f90 -- source code for IGMabs
- lmodel.dat, lmodel_igmabs.dat -- local model definition file needed by xspec.  
- compile.sh -- installation script written on bash.
- README.md -- this file

INSTALLATION

You can use the compile.sh file to install the model by doing

sh compile.sh

In the  model folder or you can setting up and using this model is as described in the xspec manual:

0) You need to have the heasoft package installed on your machine, but it must be built from source. Local models cannot be installed from a binary installation.

1) untar this directory somewhere in your user area

2) setup your headas environment (eg. 'setenv HEADAS /path/to/architecture',and 'source \$HEADAS/headas-init.csh')

3) start up xspec, and in response to the prompt type 

'initpackage igmabs lmodel_igmabs.dat <path-to-current-directory>',

where <path-to-current-directory> is the full path of the current directory. After the build is complete type 

'lmod igmabs <path-to-current-directory>'

In subsequent  sessions you don't neet to do the initpackage step again, just the lmod.

PARAMETERS

Inside of xspec, the model can be invoked by typing 'mo igmabs*pow' or variations on that. The input parameters included the elemental column densities, the turbulence velocity and redshift. 

CONTACT

This package is still being tested. Please contact me with any reports or questions.

egatuzz@mpe.mpg.de



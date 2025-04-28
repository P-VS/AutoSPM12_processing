UZB-fmr-processing-scripts

This repository contains the scripts used at UZ Brussel to standardise, optimise and automate the processing of single echo and multi-echo fMRI data. Additionally, a batch script is provided to do the segmentation in CAT12 for VBM.
The python and Matlab scripts can be used independently from each other.
The scripts comes without any warrant for the results. You may use them at your own responsibility. I only share them to help and inspire others with their own processing scripts.

The folders contains following tools:

1. Python
-convert_dcm2niix: convert dicom to nii data format using dcm2niix and organise the data for the processing (BIDS based)

2. Matlab
-my_spmbatch/AutoDCM2nii: initiating the automatic ccoonversion of dicom data into nifti data and organising them for further processing with my_spmbatch

-my_spmbatch/AutoSPMpreprocessing_fmri: initiating the automatic preprocessing of fMRI data in SPM12.

-my_spmbatch/AutoSPM1stlevel_fmri: initiating the automatic 1st level analysis of fMRI data in SPM12.

-my_spmbatch/AutoSPMpreprocessing_vbm: initiating the automatic preprocessing of VBM data (T1 anatomical high resolution scans) in CAT12.

-my_spmbatch/AutoSPMpreprocessing_asl: initiating of the automatic preprocessing of 3D PCASL data from a GE scanner (WIP)

-my_spmbatch/AutoSPMpreprocessing_aslbold: initiating the automatic preprocessing of ASLBOLD (experimental sequence on GE) (WIP)

Each analysis can be started by setting all the parameters in the beginning of the AutoSPM… script and clic ‘Run’.

3. fMRI_Class
The slides of my fMRI classes are provided. During these classes, a demo dataset is analysed manually in the same way as is done by my scripts.

Required software installed:

* Matlab
* SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
* MRIcroGL (https://www.mccauslandcenter.sc.edu/mricrogl/) for dicom to nifti conversion
* CAT12 (https://neuro-jena.github.io/cat/) for VBM
* GIFT (https://trendscenter.org/software/gift/) to use ICA-AROMA denoising

To use the Matlab script within SPM:
* ACID SPM toolbox (http://www.diffusiontools.com/index.html)
* MEHBfMRI SPM toolbox (https://github.com/P-VS/MEHBfMRI)

To use the Python script (Linux annd Mac only):
* Anaconda (https://www.anaconda.com/distribution/)
    * start Spyder in a terminal/Command Promt: type "spyder"
* Nipype (https://nipype.readthedocs.io/en/latest/):
    * in terminal type: conda install --channel conda-forge nipype
* nibabel (https://nipy.org/nibabel/):
    * In terminal type: pip install nibabel
* nilearn (https://nilearn.github.io):
    * in terminal type: pip install -U --user nilearn
* dcm2niix: conda install -c conda-forge dcm2niix

The fMRI preprocessing script for fMRI provide the following steps (In a fixed order):

IMPORTANT: !! Look at your data before starting any preprocessing. It makes no sense to lose time in trying to preprocess bad data !!

1. Automatic set the origin in the anterior comisura and align with the MNI template to improve coregistration and normalisation steps 
2. Remove dummy scans (prefix to the file name for the combination of steps 1 and 2 = e)
3. Realignment to the first echo (prefix to the file name = r)
4. TOPUP like EPI geometric distortion correction based on a phase encoding gradient polarity reversed scan (pepolar) (prefix to the file name = u)
5. For ASLBOLD: filetering between BOLD (f<0.1Hz) (prefix to the file name = f, endfix is set to _bold) and ALS (f>0.1Hz) part (endfix is set to _asl)
6. denoising (prefix to the file name = d)
    1. Bandpass filtering and detrending
    2. Extension of  motion regressors to 24 by adding the temporal derivative and the squared regressors)
    3. aCompCor (noise componenten determined on no gray or white matter voxels)
    4. ICA-AROMA for single echo fMRI / ME-ICA-AROMA foor multi-echo fMRI. A component is labeled as noise if at least 1 of following criteria is met:
        1. Fraction of the component in non brain areas > 2 * fraction of the component in brain areas (grey and white matter areas)
        2. The highest correlation between a component’s time course and the noise regressors (24 motion+aCompCor) (see Van Schuerbeek te al. The optimized combination of aCompCor and ICA-AROMA to reduce motion and physiologic noise in task fMRI data. Biomed. Physiologic. Eng. Express 2022, 8(5), doi:10.1088/2057-1976/ac63f0)
        3. Non T2*-related signal as determined with ME-ICA: Rho > 1.25 * Kappa (see Kundu et al. Differentiating BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage 2012, 60(3): 1759-1770, code based ons tedana: https://github.com/ME-ICA/tedana)
        4. High frequency content > 50%
    5. Noise regression of the 24 motion regressors, aCompCor regressors / Soft removal of the ICA noise components
7. For ASLBOLD: the S0 components as determined with ME-ICA (non ASL if Kappa > 1.25 * Rho) in step 6.3 are added to the filtered ASL component from step 5 to form the functional asl signal
8. If multi-echo fMRI: echo combination (prefix to the file name = c). Choices are
    1. Simple averaging
    2. TE weighted
    3. T2* weighted (as in tedana: https://github.com/ME-ICA/tedana) 
    4. Dynamic T2* mapping (see Heunis et al. 2021. The effects of multi-echo fMRI combination and rapid T2*-mapping on offline and real-time BOLD sensitivity. NeuroImage 238, 118244)
9. Slice time correction (also possible with HyperBand/Simultaneous multislice) (prefix to the file name = a)
10. Normalisation to the MNI template (prefix to the file name = w)
11. Smoothing (prefix to the file name = s)
12. Denoising (if not done in step 6) (prefix to the file name = d)

IMPORTANT: !! Look at your the result of the preprocessing before starting any 1st level (individual subject) or groups analysis. Doing statistical tests on wrongly preprocessed data can affect your results!

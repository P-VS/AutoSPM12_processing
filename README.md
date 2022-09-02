# UZB-fmr-processing-scripts

This repository contains some scripts used at UZ Brussel to standardise, optimise and automate the processing single echo and multi-echo fMRI data. The python and Matlab scripts can be used independently from each other.

The scripts comes without any warrant for the results. You may use them at your own responsibility. I only share them to help and inspire others with their own processing scripts.

The folders contains following tools:

1. Python

-AutoCropBET_T1w: FSL based cropping (RobustFOV), brain extraction (BET) and segmentation (FAST) of T1 weighted head scans
-AutofMRI_denoising: Denoising of fMRI timeseries. 

One can choose between noise regression, aCompCor, ICA-AROMA and combined aCompCor-ICA-AROMA (Van Schuerbeek et al. 2022, The optimized combination of aCompCor and ICA-AROMA to reduce motion and physiologic noise in task fMRI data. Biomed Phys Eng Express 1;8(5). doi: 10.1088/2057-1976/ac63f0)

-convertdcm2nii: Convert dicom to nifty format using dcm2niix (included in nipype). The resulting files are organised in BIDS format.

2. Matlab

-spm_read_vols: alternative version of spm_read_vols that loads nifty data per volume rather than per slice to fasten up data loading. This script used the scripts in the my_spmbatch/NIFTI folder (copied from https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image). To use it, you have to replace or rename the original spm_read_vols.m file in the SPM12 folder.

-my_spmbatch/AutoSPMpreprocessing: initiating the automatic preprocessing of fMRI data in SPM12.

To improve the coregistration and normalisation steps, it is advisable to set the origin and orientation of all T1w, fMRI and field map scans according to the AC-PC line.
Since some SPM preprocessing steps were found to be extremely slow with current fMRI datasets (70 sliced and 600 dynamics), I made my own faster copy of these scripts.

-my_spmbatch/AutoSPMprocessing: initiating the automatic processing of fMRI data in SPM12.

The onsets and timings of the task events/blocks are read in a events.tsv file.

If no dummy scans at the start of the fMRI scan were added, you see a decrease in the image contrast in the first few dynamics. In that case, it is advisable to delete these first dynamics from the fMRI time series. In the AutoSPMpreprocessing script a step to delete these first "dummy" scans is included and the AutoSPMprocessing script will correct the onset timings in the event.tsv file with the start of the first dynamic in the remaining fMRI data. The number of dynamics that has to be deleted is calculated based on the 'dummytime' parameter (in seconds) and the TR. If your fMRI task started after the dummy scans at the scanner, you have to set the 'dummytime' as 0.

-my_spmnatch/my_spmbatch_normfslsegmets: normalisation of the segmentation maps made in SPM (see python/AutoCropBET_T1w)

To my experience, the segmentation done in FSL FAST seems to be less dependent on the alignment (setting the origin and orientation according to AC-PC) between the T1w scan and the MNI template than SPM. However, according to Tudorascu et al. 2016 (Reproducibility and Bias in Healthy Brain Segmentation: Comparison of Two Popular Neuroimaging Platforms. Front. Neurosci., Sec. Brain Imaging Methods https://doi.org/10.3389/fnins.2016.00503) both SPM and FSL segmentation perform well, but does not give exact the same results.

--------------------------------------------------------------------------------------

Required software installed:

- Matlab
- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/)
	- make sure SPM and its subfolders are added to the Matlab Path
- MRIcroGL (https://www.mccauslandcenter.sc.edu/mricrogl/) for dicom to nifti conversion

To fasten up the preprocessing and processing of multiple fMRI data sets, Matlab's Parallel toolbox can be used.

To use the Matlab script within SPM:

- ACID SPM toolbox (http://www.diffusiontools.com/index.html)
- MEHBfMRI SPM toolbox (https://github.com/P-VS/MEHBfMRI)

To use the Python script:

- FSL
- Anaconda (https://www.anaconda.com/distribution/)
	- start Spyder in a terminal/Command Promt: type "spyder"
- DICOM Sort for sorting dicom files per series
	- from terminal: pip install DICOMsort
	- start from terminal by the command "dicomsort"
	- as program (https://dicomsort.com)
- Nipype (https://nipype.readthedocs.io/en/latest/): 
	- in terminal type: conda install --channel conda-forge nipype
- nibabel (https://nipy.org/nibabel/): 
	- In terminal type: pip install nibabel 
- nilearn (https://nilearn.github.io):
	- in terminal type: pip install -U --user nilearn
- dcm2niix: conda install -c conda-forge dcm2niix

--------------------------------------------------------------------------------------

Prior to using this script to process your data:

Convert the DICOM files into nifti using dcm2niix in MROCroGL
For the anatomical scans, set 'Crop 3D Images' on

In SPM, set the origin and orientation of all scans according to the anterior-posterior comissura

Organise the data in BIDS format
    - datpath
        -sub-##
            -ses-00# (if your experiment contains multiple session per subject)
                -anat: containes the anatomical data (3D T1)
                   Files: sub-##_T1w.nii and sub-##_T1w.json
                -func: containes the fmri data
                   Files: sub-##_task-..._bold.nii and sub-##_task-..._bold.json
                -fmap: containnes the gradient pololarity (blip-up/down) filpt data or the fieldmap scans
                   Files in case of inverted gradient polarity: sub-##_dir-pi_epi.nii and sub-##_dir-pi_epi.json
                   Files in case of fieldmap scans: (image 1 in file is amplitude, image 2 in file is phase)
                          sub-##_fmap_echo-1.nii and sub-##_fmap_echo-1.json
                          sub-##_fmap_echo-2.nii and sub-##_fmap_echo-2.json
    
IMPORTANT: !! Look at your data before starting any (pre)processing. It makes no sense to lose time in trying to process bad data !!

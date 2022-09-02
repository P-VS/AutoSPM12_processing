#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:04:03 2021

@author: dr. Peter Van Schuerbeek
"""

"""
Denoising of the fMRI data
"""


import warnings
import sys
import os

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import shutil

import nibabel as nib

from nipype import Workflow, Node, IdentityInterface

from nipype.interfaces.io import SelectFiles
from nipype.interfaces.utility import Function
from nipype.interfaces.fsl.model import MELODIC
from nipype.interfaces.fsl.utils import Smooth

"""
---------------------------------------------------------------------------------------
"""

def set_denoising_parameters():
    
    """
    Give the basic input information of your analysis
    """
    an_params = {}
    
    an_params['datpath'] = '/Volumes/LaCie/UZ_Brussel/'  #No spaties in the path name
    
    first_sub = 1
    last_sub = 5
    an_params['sublist'] = [1] #list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = list(range(first_sub,last_sub+1))
    
    an_params['nsessions'] = [1] #nsessions>0 data should be in sub-ii/ses-00j
    an_params['ME_fMRI'] = True
    
    an_params['task'] = [''] #text string that is in between task- and _bold in your fMRI nifti filename
    
    an_params['preproc_folder'] = 'preproc_func'
    
    an_params['rp_derivatives'] = True
    an_params['censoring'] = False
    an_params['CompCor'] = True
    
    an_params['do_ICA_aroma'] = True
    an_params['do_ICA_removal'] = True
    
    an_params['use_CompCor_in_AROMA'] = True
    an_params['concat_AROMA_CompCor'] = False
    
    an_params['do_noise_regression'] = False
    
    an_params['denoise_folder'] = 'preproc_func_AROMA'
    
    return an_params


"""
---------------------------------------------------------------------------------------
BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
---------------------------------------------------------------------------------------
"""
    
def compcor_dn(fmrifile,csf_file,gm_file,wm_file,t_r=2.0):
    
    import os
    
    import nibabel as nib
    
    from nilearn import masking
    
    from scipy.ndimage import binary_dilation
    from skimage.morphology import ball
    
    from nipype.algorithms.confounds import ACompCor
    
    csfim = nib.load(csf_file)
    gmim = nib.load(gm_file)
    wmim = nib.load(wm_file)
    
    gm_data = gmim.get_data()>0.05
    wm_data = wmim.get_data()
    csf_data = csfim.get_data()
    
    wm_data[wmim.get_data()<0.99]=0
    csf_data[csfim.get_data()<0.99]=0
    
    # Dilate the GM mask
    gm_data = binary_dilation(gm_data, structure=ball(2))
    
    wm_data[gm_data] = 0
    comb_data = wm_data+csf_data
    
    fmridat = nib.load(fmrifile)
    fmrimask = masking.compute_epi_mask(fmridat)
    
    boldmask_file = os.path.join(os.getcwd(),'boldmask.nii')
    
    nib.save(fmrimask,boldmask_file)
    
    wm_data[fmrimask.get_data()<0.5]=0
    csf_data[fmrimask.get_data()<0.5]=0
    comb_data[fmrimask.get_data()<0.5]=0
    
    wm_data[wm_data>0.0]=1
    csf_data[csf_data>0.0]=1
    comb_data[comb_data>0.0]=1
    
    nwmim=nib.Nifti1Image(wm_data,gmim.affine,gmim.header)
    csfmim=nib.Nifti1Image(csf_data,gmim.affine,gmim.header)
    ncombim=nib.Nifti1Image(comb_data,gmim.affine,gmim.header)
    
    wm_file = os.path.join(os.getcwd(),'compcor_wm.nii')
    csf_file = os.path.join(os.getcwd(),'compcor_csf.nii')
    comb_file = os.path.join(os.getcwd(),'compcor_csfwm.nii')
    
    nib.save(nwmim,wm_file)
    nib.save(csfmim,csf_file)
    nib.save(ncombim,comb_file)
    
    acompcor_node = ACompCor()
    acompcor_node.inputs.realigned_file = fmrifile
    acompcor_node.inputs.pre_filter='cosine'
    acompcor_node.inputs.components_file = 'acompcor_components.tsv'
    acompcor_node.inputs.header_prefix = 'acompcor_'
    acompcor_node.inputs.mask_files = [csf_file]
    acompcor_node.inputs.merge_method='none'
    #acompcor_node.inputs.variance_threshold = 0.5
    acompcor_node.inputs.num_components = 5 #'all'
    acompcor_node.inputs.repetition_time = t_r
    acompcor_node.inputs.failure_mode='NaN'
    
    acompcor_node.run()
    
    confounds_file = os.path.join(os.getcwd(),'acompcor_components.tsv')
    
    return confounds_file

"""
---------------------------------------------------------------------------------------
"""  
def mot_derivatives(rp_file,order=0,derivatives=1,fmriname=''):
    """
    Expand motion regressors upto 2nd order derivative

    motion + d(motion)/dt + d2(motion)/dt2 (linear + quadratic)
    """
    import os
    
    import numpy as np
    
    params = np.genfromtxt(rp_file)
    out_params = params
    for d in range(1, derivatives + 1):
        cparams = np.vstack((np.repeat(params[0, :][None, :], d, axis=0),params))
        out_params = np.hstack((out_params, np.diff(cparams, d, axis=0)))
    out_params2 = out_params
    for i in range(2, order + 1):
        out_params2 = np.hstack((out_params2, np.power(out_params, i)))

    rp_file = os.path.join(os.getcwd(),'confounds_'+fmriname)
    np.savetxt(rp_file, out_params2, fmt="%.10f")
    
    return rp_file

"""
---------------------------------------------------------------------------------------
"""  

def make_censor_confounds(outlier_file,tdim=1):
    import os

    import numpy as np

    outliers = np.loadtxt(outlier_file).astype(int)
    
    if outliers.shape:
        outliers_confound = np.zeros((tdim,outliers.shape[0])).astype(int)
        
        for i in range(0,outliers.shape[0]):
           outliers_confound[outliers[i],i]=1
    else:
        outliers_confound = np.zeros((tdim,1)).astype(int)
        
        outliers_confound[outliers,0]=1
       
    confounds_file = os.path.join(os.getcwd(),'censor_confounds.txt')
    
    np.savetxt(confounds_file, outliers_confound, fmt="%d")
    
    return confounds_file
                 
"""
---------------------------------------------------------------------------------------
""" 
   
def concat_confounds(file_1,file_2,fmriname=''):  
    
    import os
    
    import numpy as np
    import pandas as pd
    
    confounds_df = pd.DataFrame()
    
    if not file_1=='none':
        f1_confound_df = np.loadtxt(file_1)
        
        if f1_confound_df.size>0:
            f1_confound_df = pd.DataFrame(f1_confound_df)
    
            confounds_df = pd.concat((confounds_df, f1_confound_df), axis=1)
    
    extension = os.path.splitext(file_2)[1]
    if extension=='.txt':
        f2_confound_df = np.loadtxt(file_2)
        if f2_confound_df.size>0:
            f2_confound_df = pd.DataFrame(f2_confound_df)
    elif extension=='.tsv':
        f2_confound_df = pd.read_csv(file_2, delimiter='\t')
    
    if f2_confound_df.size>0:
        confounds_df = pd.concat((confounds_df, f2_confound_df), axis=1)
    
    confounds_file = os.path.join(os.getcwd(),'confounds_'+fmriname)
    
    np.savetxt(confounds_file, confounds_df, fmt="%.10f")
    
    return confounds_file

"""
---------------------------------------------------------------------------------------
"""

def save_preproc_files(in_file,save_dir):
    
    import os
    import shutil
    
    if len(in_file[0])>len(in_file):
        for i in range(0,len(in_file)):
            saved_file = os.path.join(save_dir,os.path.basename(in_file[i]))
            shutil.copy(in_file[i],saved_file)
    else:
        saved_file = os.path.join(save_dir,os.path.basename(in_file))
        
        if os.path.isfile(saved_file):
            os.remove(saved_file)
            
        shutil.copy(in_file,saved_file)
    
    return saved_file

"""
---------------------------------------------------------------------------------------
"""

def save_melodic_folder(in_dir,save_dir):
    
    import shutil
    
    shutil.rmtree(save_dir)
    shutil.copytree(in_dir,save_dir)
    
    return save_dir

"""
---------------------------------------------------------------------------------------
"""

def make_head_mask(anat):
    
    import os
    import nibabel as nib
    
    from nilearn import masking
    
    head_im = nib.load(anat)
    hmask = masking.compute_epi_mask(head_im)
    
    head_mask = os.path.join(os.getcwd(),'headmask.nii')
    
    nib.save(hmask,head_mask)

    return head_mask

"""
---------------------------------------------------------------------------------------
"""

"""Code largely copied from https://github.com/maartenmennes/ICA-AROMA"""


def ica_aroma(func,confounds,mel_dir,csf_file,gm_file,wm_file,head_mask,mask,t_r=2.0):
    
    def cross_correlation(a, b):
    
        import numpy as np
        
        """Cross Correlations between columns of two matrices"""
        assert a.ndim == b.ndim == 2
        _, ncols_a = a.shape
        # nb variables in columns rather than rows hence transpose
        # extract just the cross terms between cols in a and cols in b
        return np.corrcoef(a.T, b.T)[:ncols_a, ncols_a:]
    
    import os
    
    import nibabel as nib
    import random
    import numpy as np
    import pandas as pd
    
    gmim = nib.load(gm_file)
    wmim = nib.load(wm_file)
    csfim = nib.load(csf_file)
    
    headmask = nib.load(head_mask)
    
    gm_data = gmim.get_data()>0.05
    wm_data = wmim.get_data()>0.05
    csf_data = csfim.get_data()>0.05
    brain_data = gm_data+wm_data
    head_data = headmask.get_data()
    
    head_data[brain_data>0.5]=0
    head_data[csf_data>0.5]=0
    
    csf_data[csf_data<0.90]=0  
    csf_data[brain_data>0.5]=0
    
    melICim = nib.load(os.path.join(mel_dir,'melodic_IC.nii'))
    
    numIC = melICim.shape[3]
    
    """Fraction component outside GM or WM"""
    edgeFract = np.zeros(numIC)
    csfFract = np.zeros(numIC)
    for i in range(1,numIC):
        iIC=nib.load(os.path.join(mel_dir,'stats','thresh_zstat'+str(i)+'.nii'))
        
        iICdat = abs(iIC.get_data())
        
        totSum = np.sum(iICdat>0)
        nbSum = np.sum(iICdat[head_data>0]>0)
        csfSum = np.sum(iICdat[csf_data>0]>0)
        
        if totSum>0:
            edgeFract[i]=nbSum/totSum
            csfFract[i]=csfSum/totSum
        
    """High frequency content"""
    # Determine sample frequency
    Fs = 1/t_r
    
    # Determine Nyquist-frequency
    Ny = Fs/2
    
    # Load melodic_FTmix file
    FT = np.loadtxt(os.path.join(mel_dir,'melodic_FTmix'))
    
    # Determine which frequencies are associated with every row in the melodic_FTmix file  (assuming the rows range from 0Hz to Nyquist)
    f = Ny * (np.array(list(range(1, FT.shape[0] + 1)))) / (FT.shape[0])
    
    # Only include frequencies higher than 0.01Hz
    fincl = np.squeeze(np.array(np.where(f > 0.01)))
    FT = FT[fincl, :]
    f = f[fincl]
    
    # Set frequency range to [0-1]
    f_norm = (f - 0.01) / (Ny - 0.01)
    
    # For every IC; get the cumulative sum as a fraction of the total sum
    fcumsum_fract = np.cumsum(FT, axis=0) / np.sum(FT, axis=0)
    
    # Determine the index of the frequency with the fractional cumulative sum closest to 0.5
    idx_cutoff = np.argmin(np.abs(fcumsum_fract - 0.5), axis=0)
    
    # Now get the fractions associated with those indices index, these are the final feature scores
    HFC = f_norm[idx_cutoff]
    
    """Maximum robust correlation with confounds"""
    # Read melodic mix file (IC time-series) and confounds
    mix = np.loadtxt(os.path.join(mel_dir,'melodic_mix'))
    conf_model = np.loadtxt(confounds)  
    
    # Determine the maximum correlation between confounds and IC time-series
    nsplits = 1000
    nmixrows, nmixcols = mix.shape
    nrows_to_choose = int(round(0.9 * nmixrows))
    
    # Max correlations for multiple splits of the dataset (for a robust estimate)
    max_correls = np.empty((nsplits, nmixcols))
    for i in range(nsplits):
        # Select a random subset of 90% of the dataset rows (*without* replacement)
        chosen_rows = random.sample(population=range(nmixrows),k=nrows_to_choose)

        # Combined correlations between RP and IC time-series, squared and non squared
        correl_nonsquared = cross_correlation(mix[chosen_rows],
                                              conf_model[chosen_rows])
        correl_squared = cross_correlation(mix[chosen_rows]**2,
                                           conf_model[chosen_rows]**2)
        correl_both = np.hstack((correl_squared, correl_nonsquared))

        # Maximum absolute temporal correlation for every IC
        max_correls[i] = np.abs(correl_both).max(axis=1)

    # Feature score is the mean of the maximum correlation over all the random splits
    # Avoid propagating occasional nans that arise in artificial test cases
    maxRPcorr=np.nanmean(max_correls, axis=0)
    
    """ This function classifies a set of components into motion and non-motion components based on three features; 
    maximum RP correlation, high-frequency content and non-bbrain-fraction"""
    
    # Define criteria needed for classification (thresholds and hyperplane-parameters)
    thr_csf = 0.10
    thr_HFC = 0.35
    hyp = np.array([-19.9751070082159, 9.95127547670627, 24.8333160239175])
    
    # Project edge & maxRPcorr feature scores to new 1D space
    x = np.array([maxRPcorr,edgeFract])
    proj = hyp[0] + np.dot(x.T, hyp[1:])
    
    # Classify the ICs
    motionICs = np.squeeze(np.array(np.where((proj > 0) + (csfFract > thr_csf) + (HFC > thr_HFC))))
            
    funcfname = os.path.basename(func)
    funcfname = funcfname.split('.nii')[0]
        
    confounds_ICs = os.path.join(os.getcwd(),'ICA-AROMA_ICs'+funcfname+'.txt')
    
    if motionICs.size>0:
        if motionICs.ndim>0:
            np.savetxt(confounds_ICs, motionICs, fmt="%i")
        else:
            np.savetxt(confounds_ICs, np.array([int(motionICs)]), fmt="%i")
    else:
        f=open(confounds_ICs,'w')
        f.close
        
    return confounds_ICs

"""
---------------------------------------------------------------------------------------
""" 

def ica_aroma_subtraction(func,confounds_ICs,mel_dir,mask):
    
    import os
    import numpy as np
    
    from nipype.interfaces.fsl.utils import FilterRegressor
    
    import shutil
    
    funcfname = os.path.basename(func)
    funcfname = funcfname.split('.nii')[0]
    
    motionICs = np.loadtxt(confounds_ICs).astype(int)
    
    if motionICs.size>0:
            
        regfilt_node = FilterRegressor()
        regfilt_node.inputs.in_file = func
        regfilt_node.inputs.design_file = os.path.join(mel_dir,'melodic_mix')
        
        if motionICs.ndim>0:
            regfilt_node.inputs.filter_columns = list(motionICs)
        else:
            regfilt_node.inputs.filter_columns = [int(motionICs)]
            
        regfilt_node.inputs.mask = mask
        regfilt_node.inputs.output_type = 'NIFTI'
        
        regfilt_node.run()
    
        denfunc = os.path.join(os.getcwd(),funcfname+'_regfilt.nii')
        
    else:
        denfunc = os.path.join(os.getcwd(),funcfname+'_regfilt.nii')
        
        shutil.copy(func,denfunc)
        
    return denfunc

"""
---------------------------------------------------------------------------------------
"""

def make_aroma_regressors(confounds_ICs,mel_dir,fmriname,save_den_dir):
    
    import os
    import numpy as np
    import pandas as pd
    
    mix = np.loadtxt(os.path.join(mel_dir,'melodic_mix'))
    conf_model = np.loadtxt(confounds_ICs).astype(int)
    
    confounds_df = pd.DataFrame()
    
    if conf_model.size>0:
    
        aroma_regs = mix[:,conf_model]

        aroma_regs = pd.DataFrame(aroma_regs)
        confounds_df = pd.concat((confounds_df, aroma_regs), axis=1)
    
    confounds_file = os.path.join(save_den_dir,'confounds_'+fmriname)
    
    np.savetxt(confounds_file, confounds_df, fmt="%.10f")
    
    return confounds_file

"""
---------------------------------------------------------------------------------------
"""   

def confound_output_names(substring,task):
    
    confound_fname = substring+'_task-'+task+'_bold.txt'
    
    der_confounds = 'der_'+confound_fname
    censor_confounds = 'censor_'+confound_fname
    compcor_confounds = 'compcor_'+confound_fname
    aroma_confounds = 'aroma_'+confound_fname
    
    return der_confounds, censor_confounds, compcor_confounds, aroma_confounds

"""
---------------------------------------------------------------------------------------
"""

def read_json_parameter(jsonfile,parameter):

    import json
    
    with open(jsonfile,'r') as f: jsondat=json.load(f)
    f.close()
                
    value = jsondat[parameter]
    
    return value

"""
---------------------------------------------------------------------------------------
"""

def regress_confounds(confounds,funcfile,mask,t_r):
    import os
    
    import numpy as np
    import nibabel as nib
    
    from nilearn.image import clean_img
    
    confound_df = np.loadtxt(confounds)
    
    if confound_df.size>0:
        clean_func = clean_img(funcfile,confounds=confound_df,detrend=True,standardize=False,t_r=t_r, mask_img=mask)
    else:
        clean_func = clean_img(funcfile,detrend=True,standardize=False,t_r=t_r, mask_img=mask)
    
    res_func =  os.path.join(os.getcwd(),'d'+os.path.basename(funcfile))
    
    nib.save(clean_func,res_func)
    
    return res_func

"""
---------------------------------------------------------------------------------------
"""

def sum_mean_denoiefunc(funcfile,dfuncfile):
    import os
    
    import nibabel as nib
    
    from nilearn.image import mean_img
    
    clean_func = nib.load(dfuncfile)
    mean_func = mean_img(funcfile)
    
    cleandat = clean_func.get_data()
    meandat = mean_func.get_data()
    
    for i in range(0,cleandat.shape[3]):
        cleandat[:,:,:,i]=cleandat[:,:,:,i]+meandat[:,:,:]
        
    resultim = nib.Nifti1Image(cleandat,clean_func.affine,clean_func.header)
    
    res_func =  os.path.join(os.getcwd(),'m'+os.path.basename(dfuncfile))
    
    nib.save(resultim,res_func)
    
    return res_func
                    
"""
---------------------------------------------------------------------------------------
"""

def main():
                   
    """
    Give the basic input information of your data
    """
    an_params = set_denoising_parameters()
    
    datpath = an_params['datpath']

    sublist = an_params['sublist']
    
    nsessions = an_params['nsessions']
    
    task = an_params['task']
    
    preproc_folder = an_params['preproc_folder']
    
    rp_derivatives = an_params['rp_derivatives']
    censoring = an_params['censoring']
    CompCor = an_params['CompCor']
    
    do_ICA_aroma = an_params['do_ICA_aroma']
    do_ICA_removal = an_params['do_ICA_removal']
    
    use_CompCor_in_AROMA = an_params['use_CompCor_in_AROMA']
    concat_AROMA_CompCor = an_params['concat_AROMA_CompCor']
    
    if use_CompCor_in_AROMA or concat_AROMA_CompCor:
        CompCor = True
        
    if use_CompCor_in_AROMA:
        concat_AROMA_CompCor = False
    
    do_noise_regression = an_params['do_noise_regression']
    
    if do_ICA_aroma:
        rp_derivatives = True
    
    
    """
    ---------------------------------------------------------------------------------------
    """
  
    print('Start preprocessing of the data')
    
    sesstring = list()
    
    do_melodic = False
    
    ferror = 0
        
    for j in nsessions:
        sesstring.append('ses-00'+str(j))
        
        for k in range(0,len(task)):
            
            substringslist = list()
            
            for i in sublist:
                if i<10:
                    substring = 'sub-0'+str(i)
                else:
                    substring = 'sub-'+str(i)
                    
                substringslist.append(substring)
            
                subpath = os.path.join(datpath,substring)
                
                subpath = os.path.join(subpath,'ses-00'+str(j))
                
                subfmridir = os.path.join(subpath,preproc_folder)
                subanadir = os.path.join(subpath,'preproc_anat')
        
                if not an_params['ME_fMRI']:
                    jsonf = os.path.join(subpath,'func',substring+'_task-'+task[k]+'_bold.json')
                    
                    subfmridat = os.path.join(subfmridir,'wraut'+substring+'_task-'+task[k]+'_bold.nii')
                    ssubfmridat = os.path.join(subfmridir,'swraut'+substring+'_task-'+task[k]+'_bold.nii')
                    subrpdat = os.path.join(subfmridir,'rp_t'+substring+'_task-'+task[k]+'_bold.txt')
                    
                    outldat = os.path.join(subfmridir,'art.ut'+substring+'_task-'+task[k]+'_bold_outliers.txt')
                    maskdat = os.path.join(subfmridir,'mask_swraut'+substring+'_task-'+task[k]+'_bold.nii')
                else:
                    jsonf = os.path.join(subpath,'func',substring+'_task-'+task[k]+'_bold_e1.json')
                    subfmridat = os.path.join(subfmridir,'wraucet'+substring+'_task-'+task[k]+'_bold.nii')
                    ssubfmridat = os.path.join(subfmridir,'swraucet'+substring+'_task-'+task[k]+'_bold.nii')
                    subrpdat = os.path.join(subfmridir,'rp_cet'+substring+'_task-'+task[k]+'_bold.txt')
                    
                    outldat = os.path.join(subfmridir,'art.ucet'+substring+'_task-'+task[k]+'_bold_outliers.txt')
                    maskdat = os.path.join(subfmridir,'mask_swraucet'+substring+'_task-'+task[k]+'_bold.nii')
   
                subanadat = os.path.join(subanadir,'wr'+substring+'_T1w_ROI.nii')
                subcbfdat = os.path.join(subanadir,'wr'+substring+'_T1w_ROI_brain_pve_0.nii')
                subgmdat = os.path.join(subanadir,'wr'+substring+'_T1w_ROI_brain_pve_1.nii')
                subwmdat = os.path.join(subanadir,'wr'+substring+'_T1w_ROI_brain_pve_2.nii')
                
                """
                Check the existence of all files
                """
            
                if not os.path.isdir(subpath):
                    print('Directory '+subpath+' not found.')
                    ferror = ferror+1
                if not os.path.isdir(subfmridir):
                    print('Directory '+subfmridir+' not found.')
                    ferror = ferror+1
                if not os.path.isdir(subanadir):
                    print('Directory '+subanadir+' not found.')
                    ferror = ferror+1
            
                if not os.path.isfile(subanadat):
                    print('File '+subanadat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(subcbfdat):
                    print('File '+subcbfdat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(subgmdat):
                    print('File '+subgmdat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(subwmdat):
                    print('File '+subwmdat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(subfmridat):
                    print('File '+subfmridat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(ssubfmridat):
                    print('File '+ssubfmridat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(subrpdat):
                    print('File '+subrpdat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(jsonf):
                    print('File '+jsonf+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(outldat):
                    print('File '+outldat+' not found.')
                    ferror = ferror+1
                if not os.path.isfile(maskdat):
                    print('File '+maskdat+' not found.')
                    ferror = ferror+1
                
                if do_ICA_aroma:
                    save_melodic_dir = os.path.join(subpath,'melodic_'+'wraut'+substring+'_task-'+task[k]+'_bold')
                    
                    if not os.path.isdir(save_melodic_dir) \
                        or not os.path.isfile(os.path.join(save_melodic_dir,'melodic_IC.nii')):
                            if not os.path.isdir(save_melodic_dir): os.mkdir(save_melodic_dir)
                            do_melodic = True

                save_den_dir = os.path.join(subpath,an_params['denoise_folder'])
                
                if not os.path.isdir(save_den_dir): os.mkdir(save_den_dir) 
        
    """
    If all data exist, start processing the data
    """
    if ferror == 0:
            
        templates = {}
        templates['anat'] = os.path.join(datpath,'{substring}','{sesstring}','preproc_anat','wr'+'{substring}'+'_T1w_ROI.nii')
        templates['csf'] = os.path.join(datpath,'{substring}','{sesstring}','preproc_anat','wr'+'{substring}'+'_T1w_ROI_brain_pve_0.nii')
        templates['gm'] = os.path.join(datpath,'{substring}','{sesstring}','preproc_anat','wr'+'{substring}'+'_T1w_ROI_brain_pve_1.nii')
        templates['wm'] = os.path.join(datpath,'{substring}','{sesstring}','preproc_anat','wr'+'{substring}'+'_T1w_ROI_brain_pve_2.nii')
        
        if not an_params['ME_fMRI']:
            templates['json'] = os.path.join(datpath,'{substring}','{sesstring}','func','{substring}'+'_task-'+'{task}'+'_bold.json')
            
            templates['func'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'wraut'+'{substring}'+'_task-'+'{task}'+'_bold.nii')
            #templates['sfunc'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'swraut'+'{substring}'+'_task-'+'{task}'+'_bold.nii')
            templates['rp_file'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'rp_t'+'{substring}'+'_task-'+'{task}'+'_bold.txt')
            
            templates['outliers'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'art.ut'+'{substring}'+'_task-'+'{task}'+'_bold_outliers.txt')
            templates['mask'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'mask_swraut'+'{substring}'+'_task-'+'{task}'+'_bold.nii')
     
        else:
            templates['json'] = os.path.join(datpath,'{substring}','{sesstring}','func','{substring}'+'_task-'+'{task}'+'_bold_e1.json')
            
            templates['func'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'wraucet'+'{substring}'+'_task-'+'{task}'+'_bold.nii')
            #templates['sfunc'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'swraucet'+'{substring}'+'_task-'+'{task}'+'_bold.nii')
            templates['rp_file'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'rp_cet'+'{substring}'+'_task-'+'{task}'+'_bold.txt')
     
            templates['outliers'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'art.ucet'+'{substring}'+'_task-'+'{task}'+'_bold_outliers.txt')
            templates['mask'] = os.path.join(datpath,'{substring}','{sesstring}',preproc_folder,'mask_swraucet'+'{substring}'+'_task-'+'{task}'+'_bold.nii')
        
        if do_ICA_aroma:
            templates['melodic_dir'] = os.path.join(datpath,'{substring}','{sesstring}','melodic_'+'wraut'+'{substring}'+'_task-'+'{task}'+'_bold')
            
        templates['save_den_dir'] = os.path.join(datpath,'{substring}','{sesstring}',an_params['denoise_folder'])
                        
        """
        Create a preprocessing workflow and select input files
        """
        infosource = Node(IdentityInterface(fields=['substring','sesstring','task']),name='infosource')
    
        infosource.iterables = [('substring', substringslist),('sesstring',sesstring),('task',task)]
    
        preproc = Workflow(base_dir=datpath,name='denoise')        

        selectfiles = Node(SelectFiles(templates,base_directory=datpath),name="selectfiles")
        
        preproc.connect(infosource, 'substring', selectfiles, 'substring')
        preproc.connect(infosource, 'sesstring', selectfiles, 'sesstring')
        preproc.connect(infosource, 'task', selectfiles, 'task')
        
        conf_names = Node(interface=Function(input_names=['substring','task'],
                                             output_names=['der_confounds','censor_confounds','compcor_confounds','aroma_confounds'],
                                             function=confound_output_names),name='confound_names')
        
        preproc.connect(infosource, 'substring', conf_names, 'substring')
        preproc.connect(infosource, 'task', conf_names, 'task')
        
        read_tr = Node(interface=Function(input_names=['jsonfile','parameter'],
                                          output_names=['value'],
                                          function=read_json_parameter),name='read_tr')
        
        read_tr.inputs.parameter = "RepetitionTime"
        
        preproc.connect([(selectfiles,read_tr,[('json','jsonfile')])])
                  
        """
        Motion derivative
        """
        if rp_derivatives:
            
            exrp_node = Node(interface=Function(input_names=['rp_file','order','derivatives','fmriname'],
                                                output_names=['rp_file'],
                                                function=mot_derivatives),name='rp_derivatives')
            
            exrp_node.inputs.order = 2 #quadratic expansion = 2
            exrp_node.inputs.derivatives = 1
            
            save_exrp_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                     output_names=['saved_file'],
                                                     function=save_preproc_files),name='save_exrp')
            
            preproc.connect([(selectfiles,exrp_node,[('rp_file','rp_file')]),
                             (conf_names,exrp_node,[('der_confounds','fmriname')]),
                             (exrp_node,save_exrp_node,[('rp_file','in_file')]),
                             (selectfiles,save_exrp_node,[('save_den_dir','save_dir')])
                             ])
                    
        """
        Censoring
        """
        if censoring:
                
            censconf_node = Node(interface=Function(input_names=['outlier_file','tdim'],
                                                    output_names=['confounds_file'],
                                                    function=make_censor_confounds),name='censor_confounds')
            
            fdat=nib.load(subfmridat)
            censconf_node.inputs.tdim = fdat.shape[3]
            
            concat_censor_node = Node(interface=Function(input_names=['file_1','file_2','fmriname'],
                                                             output_names=['confounds_file'],
                                                             function=concat_confounds),name='concat_censor')
            
            save_cens_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                     output_names=['saved_file'],
                                                     function=save_preproc_files),name='save_cens')
            
            if rp_derivatives:
                preproc.connect([(selectfiles,censconf_node,[('outliers','outlier_file')]),
                                 (conf_names,concat_censor_node,[('censor_confounds','fmriname')]),
                                 (exrp_node,concat_censor_node,[('rp_file','file_1')]),
                                 (censconf_node,concat_censor_node,[('confounds_file','file_2')]),
                                 (concat_censor_node,save_cens_node,[('confounds_file','in_file')]),
                                 (selectfiles,save_cens_node,[('save_den_dir','save_dir')])
                                ])
            else:
                preproc.connect([(selectfiles,censconf_node,[('outliers','outlier_file')]),
                                 (conf_names,concat_censor_node,[('censor_confounds','fmriname')]),
                                 (selectfiles,concat_censor_node,[('rp_file','file_1')]),
                                 (censconf_node,concat_censor_node,[('confounds_file','file_2')]),
                                 (concat_censor_node,save_cens_node,[('confounds_file','in_file')]),
                                 (selectfiles,save_cens_node,[('save_den_dir','save_dir')])
                                ])
                        
        """
        aCompCor
        """
        if CompCor:
            compcor_node=Node(interface=Function(input_names=['fmrifile','csf_file','gm_file','wm_file'],
                                                 output_names=['confounds_file'],
                                                 function=compcor_dn),name='compcor')
            
            preproc.connect([(selectfiles,compcor_node,[('func','fmrifile'),
                                                        ('csf','csf_file'),
                                                        ('gm','gm_file'),
                                                        ('wm','wm_file'),
                                                        ]),
                             (read_tr,compcor_node,[('value','t_r')])
                            ])
            
            save_compcor_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                        output_names=['saved_file'],
                                                        function=save_preproc_files),name='save_compcor')
        
            concat_compcor_node = Node(interface=Function(input_names=['file_1','file_2','fmriname'],
                                                      output_names=['confounds_file'],
                                                      function=concat_confounds),name='concat_rp-ccor')
            
            if censoring:                
                preproc.connect([(concat_censor_node,concat_compcor_node,[('confounds_file','file_1')]),
                                 (compcor_node,concat_compcor_node,[('confounds_file','file_2')]),
                                 (conf_names,concat_compcor_node,[('compcor_confounds','fmriname')]),
                                 (concat_compcor_node,save_compcor_node,[('confounds_file','in_file')]),
                                 (selectfiles,save_compcor_node,[('save_den_dir','save_dir')])
                                 ])
            elif rp_derivatives:                
                preproc.connect([(exrp_node,concat_compcor_node,[('rp_file','file_1')]),
                                 (compcor_node,concat_compcor_node,[('confounds_file','file_2')]),
                                 (conf_names,concat_compcor_node,[('compcor_confounds','fmriname')]),
                                 (concat_compcor_node,save_compcor_node,[('confounds_file','in_file')]),
                                 (selectfiles,save_compcor_node,[('save_den_dir','save_dir')])
                                 ])
            else:                
                preproc.connect([(selectfiles,concat_compcor_node,[('rp_file','file_1')]),
                                 (compcor_node,concat_compcor_node,[('confounds_file','file_2')]),
                                 (conf_names,concat_compcor_node,[('compcor_confounds','fmriname')]),
                                 (concat_compcor_node,save_compcor_node,[('confounds_file','in_file')]),
                                 (selectfiles,save_compcor_node,[('save_den_dir','save_dir')])
                                 ])
            
        """
        ICA-AROMA
        """
        if do_ICA_aroma:
            head_mask_node = Node(interface=Function(input_names=['anat'],
                                                     output_names=['head_mask'],
                                                     function=make_head_mask),name='head_mask')
            
            preproc.connect([(selectfiles,head_mask_node,[('anat','anat')])])
            
            if do_melodic:
                
                melodic_node = Node(MELODIC(output_type='NIFTI',out_all=True),name='melodic')
                
                save_melodic_node = Node(interface=Function(input_names=['in_dir','save_dir'],
                                                            output_names=['save_dir'],
                                                            function=save_melodic_folder),name='save_melodic')
        
                preproc.connect([(selectfiles,melodic_node,[('func','in_files')]),
                                 (head_mask_node,melodic_node,[('head_mask','mask')]),
                                 (read_tr,melodic_node,[('value','tr_sec')]),
                                 (melodic_node,save_melodic_node,[('out_dir','in_dir')]),
                                 (selectfiles,save_melodic_node,[('melodic_dir','save_dir')])
                                 ])
            
            
            ica_aroma_node = Node(interface=Function(input_names=['func','confounds','mel_dir','csf_file','gm_file','wm_file','head_mask','mask','t_r'],
                                                     output_names=['confounds_ICs'],
                                                     function=ica_aroma),name='ica_aroma')
            
            save_confIC_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                       output_names=['saved_file'],
                                                       function=save_preproc_files),name='save_confIC')
            
            if use_CompCor_in_AROMA:
                preproc.connect([(concat_compcor_node,ica_aroma_node,[('confounds_file','confounds')])])             
            else:
                preproc.connect([(exrp_node,ica_aroma_node,[('rp_file','confounds')])])
            
            if do_melodic:
                preproc.connect([(melodic_node,ica_aroma_node,[('out_dir','mel_dir')])])
            else:
                preproc.connect([(selectfiles,ica_aroma_node,[('melodic_dir','mel_dir')])])
                
            preproc.connect([(selectfiles,ica_aroma_node,[('func','func'),#('sfunc','func'),
                                                          ('mask','mask'),
                                                          ('csf','csf_file'),
                                                          ('gm','gm_file'),
                                                          ('wm','wm_file')
                                                         ]),
                             (read_tr,ica_aroma_node,[('value','t_r')]),
                             (head_mask_node,ica_aroma_node,[('head_mask','head_mask')]),
                             (ica_aroma_node,save_confIC_node,[('confounds_ICs','in_file')]),
                             (selectfiles,save_confIC_node,[('save_den_dir','save_dir')])
                             ])   
            
            if do_noise_regression:
                icareg_node=Node(interface=Function(input_names=['confounds_ICs','mel_dir','fmriname','save_den_dir'],
                                                    output_names=['confounds_file'],
                                                    function=make_aroma_regressors),name='make_aroma_regressors')
                
                preproc.connect([(ica_aroma_node,icareg_node,[('confounds_ICs','confounds_ICs')]),
                                 (conf_names,icareg_node,[('aroma_confounds','fmriname')]),
                                 (selectfiles,icareg_node,[('save_den_dir','save_den_dir')])
                                 ])
                
                if do_melodic:
                    preproc.connect([(melodic_node,icareg_node,[('out_dir','mel_dir')])])
                else:
                    preproc.connect([(selectfiles,icareg_node,[('melodic_dir','mel_dir')])])
                    
                if concat_AROMA_CompCor:
                    concat_aroma_node = Node(interface=Function(input_names=['file_1','file_2','fmriname'],
                                                                output_names=['confounds_file'],
                                                                function=concat_confounds),name='concat_aroma-ccor')
                    
                    save_aromacc_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                                output_names=['saved_file'],
                                                                function=save_preproc_files),name='save_aroma-cc')
                    
                    preproc.connect([(icareg_node,concat_aroma_node,[('confounds_file','file_1')]),
                                     (compcor_node,concat_aroma_node,[('confounds_file','file_2')]),
                                     (conf_names,concat_aroma_node,[('aroma_confounds','fmriname')]),
                                     (concat_aroma_node,save_aromacc_node,[('confounds_file','in_file')]),
                                     (selectfiles,save_aromacc_node,[('save_den_dir','save_dir')])
                                     ])
                        
                
        """
        ICA-AROMA noise component subtraction
        """
        if do_ICA_removal:
            icasub_node=Node(interface=Function(input_names=['func','confounds_ICs','mel_dir','mask'],
                                                output_names=['denfunc'],
                                                function=ica_aroma_subtraction),name='ica_removal')
            
            save_icasub_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                       output_names=['saved_file'],
                                                       function=save_preproc_files),name='save_ica_aroma')
            
            preproc.connect([(selectfiles,icasub_node,[('func','func'),#('sfunc','func'),
                                                       ('mask','mask')]),
                             (ica_aroma_node,icasub_node,[('confounds_ICs','confounds_ICs')]),
                             (icasub_node,save_icasub_node,[('denfunc','in_file')]),
                             (selectfiles,save_icasub_node,[('save_den_dir','save_dir')])
                             ])
            
            if do_melodic:
                preproc.connect([(melodic_node,icasub_node,[('out_dir','mel_dir')])])
            else:
                preproc.connect([(selectfiles,icasub_node,[('melodic_dir','mel_dir')])])
                        
                   
        """
        Noise regression
        """
        if do_noise_regression:
            regconf_node=Node(interface=Function(input_names=['confounds','funcfile','mask','t_r'],
                                                 output_names=['res_func'],
                                                 function=regress_confounds),name='regconf')
            
            """save_reconf_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                       output_names=['saved_file'],
                                                       function=save_preproc_files),name='save_regconf')"""
            
            sum_meandn_node=Node(interface=Function(input_names=['funcfile','dfuncfile'],
                                                    output_names=['res_func'],
                                                    function=sum_mean_denoiefunc),name='summeandf')
            
            save_meandn_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                       output_names=['saved_file'],
                                                       function=save_preproc_files),name='save_meandn')
            
            preproc.connect([(selectfiles,regconf_node,[('func','funcfile')]),#[('sfunc','funcfile')]),
                             (read_tr,regconf_node,[('value','t_r')])
                             ])
            
            if concat_AROMA_CompCor:
                preproc.connect([(selectfiles,regconf_node,[('mask','mask')]),
                                 (concat_aroma_node,regconf_node,[('confounds_file','confounds')]),
                                 (selectfiles,sum_meandn_node,[('sfunc','funcfile')]),
                                 (regconf_node,sum_meandn_node,[('res_func','dfuncfile')]),
                                 (sum_meandn_node,save_meandn_node,[('res_func','in_file')]),
                                 (selectfiles,save_meandn_node,[('save_den_dir','save_dir')])
                                 ])
            elif do_ICA_aroma:
                preproc.connect([(selectfiles,regconf_node,[('mask','mask')]),
                                 (icareg_node,regconf_node,[('confounds_file','confounds')]),
                                 (selectfiles,sum_meandn_node,[('func','funcfile')]),#(selectfiles,sum_meandn_node,[('sfunc','funcfile')]),
                                 (regconf_node,sum_meandn_node,[('res_func','dfuncfile')]),
                                 (sum_meandn_node,save_meandn_node,[('res_func','in_file')]),
                                 (selectfiles,save_meandn_node,[('save_den_dir','save_dir')])
                                 ])
            elif CompCor:
                preproc.connect([(selectfiles,regconf_node,[('mask','mask')]),
                                 (concat_compcor_node,regconf_node,[('confounds_file','confounds')]),
                                 (selectfiles,sum_meandn_node,[('func','funcfile')]),#(selectfiles,sum_meandn_node,[('sfunc','funcfile')]),
                                 (regconf_node,sum_meandn_node,[('res_func','dfuncfile')]),
                                 (sum_meandn_node,save_meandn_node,[('res_func','in_file')]),
                                 (selectfiles,save_meandn_node,[('save_den_dir','save_dir')])
                                 ])
            
            elif censoring:
                preproc.connect([(selectfiles,regconf_node,[('mask','mask')]),
                                 (concat_censor_node,regconf_node,[('confounds_file','confounds')]),
                                 (selectfiles,sum_meandn_node,[('func','funcfile')]),#(selectfiles,sum_meandn_node,[('sfunc','funcfile')]),
                                 (regconf_node,sum_meandn_node,[('res_func','dfuncfile')]),
                                 (sum_meandn_node,save_meandn_node,[('res_func','in_file')]),
                                 (selectfiles,save_meandn_node,[('save_den_dir','save_dir')])
                                 ])
            elif rp_derivatives:
                preproc.connect([(selectfiles,regconf_node,[('mask','mask')]),
                                 (exrp_node,regconf_node,[('rp_file','confounds')]),
                                 (selectfiles,sum_meandn_node,[('func','funcfile')]),#(selectfiles,sum_meandn_node,[('sfunc','funcfile')]),
                                 (regconf_node,sum_meandn_node,[('res_func','dfuncfile')]),
                                 (sum_meandn_node,save_meandn_node,[('res_func','in_file')]),
                                 (selectfiles,save_meandn_node,[('save_den_dir','save_dir')])
                                 ])
            else:
                preproc.connect([(selectfiles,regconf_node,[('mask','mask')]),
                                 (selectfiles,regconf_node,[('rp_file','confounds')]),
                                 (selectfiles,sum_meandn_node,[('func','funcfile')]),#(selectfiles,sum_meandn_node,[('sfunc','funcfile')]),
                                 (regconf_node,sum_meandn_node,[('res_func','dfuncfile')]),
                                 (sum_meandn_node,save_meandn_node,[('res_func','in_file')]),
                                 (selectfiles,save_meandn_node,[('save_den_dir','save_dir')])
                                 ])
                
        """
        FSL Smooth
        """
        
        smooth_node = Node(Smooth(fwhm=6, output_type='NIFTI'),name='smooth')
        
        save_smooth_node = Node(interface=Function(input_names=['in_file','save_dir'],
                                                   output_names=['saved_file'],
                                                   function=save_preproc_files),name='save_smooth')
        
        if do_ICA_removal:
            preproc.connect([(icasub_node,smooth_node,[('denfunc','in_file')]),
                             (smooth_node,save_smooth_node,[('smoothed_file','in_file')]),
                             (selectfiles,save_smooth_node,[('save_den_dir','save_dir')])
                             ])
        else:   
            preproc.connect([(sum_meandn_node,smooth_node,[('res_func','in_file')]),
                             (smooth_node,save_smooth_node,[('smoothed_file','in_file')]),
                             (selectfiles,save_smooth_node,[('save_den_dir','save_dir')])
                             ])
                            
        """
        Run the workflow
        """
        
        print('')
        print('Start denoising')
        print('')
        preproc.run(plugin='MultiProc')
        #preproc.run()
        print('Done denoising')
        print('')
        
        """
        #Clean up intermediate saved data
        """
        shutil.rmtree(os.path.join(datpath,'denoise'), ignore_errors=True)
                
if __name__ == '__main__':
    main()

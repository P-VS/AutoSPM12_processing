#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 13:14:04 2021

@author: dr. Peter Van Schuerbeek
"""

"""
Converting DICOM into .nii data

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
                -perf: ASL data based on GE PCASL sequence which gives only the m0, deltam and CBF images
                    Files: sub-##_asl_m0scan.nii
                    Files: sub-##_asl_deltam.nii
                    Files: sub-##_asl_cbf.nii
"""

import warnings
import sys
import os

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import shutil
import json

from nipype import Workflow, Node, IdentityInterface

from nipype.interfaces.io import SelectFiles
from nipype.interfaces.utility import Function
from nipype.interfaces.dcm2nii import Dcm2niix

def set_preprocessing_parameters():
    
    """
    Give the basic input information of your data
    """
    pp_params = {}
    
    pp_params['datpath'] = '/Volumes/LaCie/UZ_Brussel/rTMS-fMRI_Interleaved/Data'  #No spaties in the path name
    
    first_sub = 1
    last_sub = 1
    pp_params['sublist'] = [3] #list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = list(range(first_sub,last_sub+1))
    
    #Add per sequence to convert an extra ssequence object to the mri_data structure as (folder,seqtype,name,[session])
    #folder: substructure starting from sub-##
    #seqtype: (anat, func, fmap, dti, asl) will be used to organise the data in folders
    #name: name of the sequence used in sub-##_... .nii to add series number to the name use %s, to add echo number use %e
    #session: scan session
    
    pp_params['mri_data'] = []
    
    scan_1 = ('dicom/T1w','anat','T1w',[1])
    #scan_2 = ('dicom/SEFMRI_PI','fmap','dir-pi_epi',[1])
    #scan_3 = ('DCM/ME-fMRI_PI','fmap','dir-pi_epi',[1])
    #scan_4 = ('DCM/SE-fMRI_EFT','func','task-SE-EFT_bold',[1])
    #scan_5 = ('dicom/SEFMRI_EMOFACES','func','task-SE-EmoFaces_bold',[1])
    #scan_6 = ('DCM/ME-fMRI_EFT','func','task-ME-EFT_bold_e%e',[1])
    #scan_7 = ('dicom/voor/fmri_task','func','task-_%s',[1])
    #scan_8 = ('dicom/voor/asl','perf','asl-1',[1])
    #scan_9 = ('dicom/tijdens/','perf','asl-%s',[1])
    #scan_10 = ('dicom/na/asl','perf','asl-5',[2])
    
    pp_params['mri_data'].append([scan_1])
    
    return pp_params


"""
BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
---------------------------------------------------------------------------------------
"""

def load_d2n_params(subdir,folder,seqtype,name,session):
    
    import os
    
    split_path = subdir.split(os.sep)
    substring = split_path[len(split_path)-1]
    
    sesdir = 'ses-00'+str(session)
    
    source_dir = os.path.join(subdir,folder)
    output_dir = os.path.join(subdir,sesdir,seqtype)
    
    out_filename = substring+'_'+name
    
    crop = False
    ignore_deriv = False
    
    if 'anat' in seqtype:
        crop = True
        ignore_deriv = True
        
    if 'fmap' in seqtype:
        if 'echo-' in name:
            out_filename = out_filename+'%s'
    
    return source_dir, crop, ignore_deriv, out_filename, output_dir

"""
---------------------------------------------------------------------------------------
"""

def main():
    
    pp_params = set_preprocessing_parameters()
    
    datpath = pp_params['datpath']

    sublist = pp_params['sublist']
    
    mri_data = pp_params['mri_data']
    
    substringslist = list()
    
    for i in sublist:

        if i<10:
            substring = 'sub-0'+str(i)
        else:
            substring = 'sub-'+str(i)
            
        substringslist.append(substring)
        
        for k in range(0,len(mri_data[0])):
            output_dir = os.path.join(datpath,substring,'ses-00'+str(mri_data[0][k][3][0]))
            if not os.path.isdir(output_dir): os.mkdir(output_dir) 
            
            output_dir = os.path.join(datpath,substring,'ses-00'+str(mri_data[0][k][3][0]),mri_data[0][k][1])
            if not os.path.isdir(output_dir): os.mkdir(output_dir) 
            
    
    templates = {}
    
    templates['subdir'] = os.path.join(datpath,'{substring}')
    
    """
    Create a preprocessing workflow and select input files
    """
    print('Make workflow step: initiate')
    
    infosource = Node(IdentityInterface(fields=['substring']),name='infosource')

    infosource.iterables = [('substring', substringslist)]

    dcm2niiwf = Workflow(base_dir=datpath,name='dcm2niiwf')
    
    selectfiles = Node(SelectFiles(templates,base_directory=datpath),name="selectfiles")
    
    dcm2niiwf.connect(infosource, 'substring', selectfiles, 'substring')
    
    for k in range(0,len(mri_data[0])):
        load_par_node = Node(interface=Function(input_names=['subdir','folder','seqtype','name','session'],
                                                output_names=['source_dir','crop','ignore_deriv','out_filename','output_dir'],
                                                function=load_d2n_params),name='load_par'+str(k))
        
        load_par_node.inputs.folder = mri_data[0][k][0]
        load_par_node.inputs.seqtype = mri_data[0][k][1]
        load_par_node.inputs.name = mri_data[0][k][2]
        load_par_node.inputs.session = mri_data[0][k][3][0]
        
        
        dcm2niiwf.connect([(selectfiles,load_par_node,[('subdir','subdir')])])
        
        d2nii_node = Node(Dcm2niix(compress='n', anon_bids=True, bids_format=True), name='d2nii'+str(k))
        
        dcm2niiwf.connect([(load_par_node,d2nii_node,[('source_dir','source_dir'),
                                                      ('crop','crop'),
                                                      ('ignore_deriv','ignore_deriv'),
                                                      ('out_filename','out_filename'),
                                                      ('output_dir','output_dir')
                                                      ])
                           ])
    
    """
    Run the workflow
    """

    print('')
    print('Start dcm2nii sub '+str(i))
    print('')
    dcm2niiwf.run()
    #dcm2niiwf.run(plugin='MultiProc')
    print('Done dcm2nii sub '+str(i))
    print('')

    shutil.rmtree(os.path.join(datpath,'dcm2niiwf'), ignore_errors=True)
    
    for k in range(0,len(mri_data[0])):
        if 'fmap' in mri_data[0][k][1]:
            if 'echo-' in mri_data[0][k][2]:
                for i in sublist:
                    if i<10:
                        substring = 'sub-0'+str(i)
                    else:
                        substring = 'sub-'+str(i)
                        
                    fmap_folder = os.path.join(datpath,substring,'ses-00'+str(mri_data[0][k][3][0]),'fmap')
                    
                    dirlist = sorted([fn for fn in os.listdir(fmap_folder) if fn.endswith('.nii')
                                                                              and not fn.startswith('.') 
                                                                              and os.path.isfile(os.path.join(fmap_folder,fn))
                                      ])
                    
                    for file in dirlist:
                        if not '_ph' in file:
                            sfile = file.split('.nii')
                            jfile = sfile[0]+'.json'
                            nfile = sfile[0]+'_am.nii'
                            njfile = sfile[0]+'_am.json'
                            
                            os.rename(os.path.join(fmap_folder,file),os.path.join(fmap_folder,nfile))
                            os.rename(os.path.join(fmap_folder,jfile),os.path.join(fmap_folder,njfile))
                            
        if 'perf' in mri_data[0][k][1]:
            for i in sublist:
                if i<10:
                    substring = 'sub-0'+str(i)
                else:
                    substring = 'sub-'+str(i)
                    
                perf_folder = os.path.join(datpath,substring,'ses-00'+str(mri_data[0][k][3][0]),'perf')
                
                dirlist = sorted([fn for fn in os.listdir(perf_folder) if fn.endswith('.nii')
                                                                          and not fn.startswith('.') 
                                                                          and os.path.isfile(os.path.join(perf_folder,fn))
                                  ])
                
                for file in dirlist:
                    if not '.json' in file:
                        sfile = file.split('.nii')
                        jfile = sfile[0]+'.json'
                        
                        with open(os.path.join(perf_folder,jfile),'r') as f: jsondat=json.load(f)
                        f.close()
                                    
                        ImType = jsondat["ImageType"]
                        
                        if 'ORIGINAL' in ImType[0]:
                            ftype = 'm0scan'
                        elif 'DERIVED' in ImType[0]:
                            if 'PRIMARY' in ImType[1]:
                                ftype = 'deltam'
                            elif 'SECONDARY' in ImType[1]:
                                ftype = 'cbf'
                         
                        if '%s' in mri_data[0][k][2]:
                            snum = '-'+str(jsondat["SeriesNumber"])+'_'
                        else:
                            snum = '_'    
                                
                        sfname = file.split('asl')
                        nfname = sfname[0]+'asl'+snum+ftype
                        
                        os.rename(os.path.join(perf_folder,file),os.path.join(perf_folder,nfname+'.nii'))
                        os.rename(os.path.join(perf_folder,jfile),os.path.join(perf_folder,nfname+'.json'))
            
                    
if __name__ == '__main__':
    main()
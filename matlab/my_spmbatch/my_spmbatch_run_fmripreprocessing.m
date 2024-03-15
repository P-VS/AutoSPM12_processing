function out = my_spmbatch_run_fmripreprocessing(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

%try
    %% preprocess anatomical scans
    if params.preprocess_anatomical
        [delfiles,keepfiles] = my_spmbatch_preprocess_anat(sub,ses,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params.save_intermediate_results,'anat','preproc_anat');
    end
    
    %% preprocess functional scans
    if params.preprocess_functional
        [delfiles,keepfiles] = my_spmbatch_functional(sub,ses,run,task,datpath,params);
        
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params.save_intermediate_results,'func',params.save_folder);  
    end
    
    %% denoise functional scans
    if params.do_denoising
        [delfiles,keepfiles] = my_spmbatch_denoise(sub,ses,run,task,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params.save_intermediate_results,params.save_folder,params.save_folder);
    end
%catch e
%    fprintf('\nPP_Error\n');
%    fprintf('\nThe error was: \n%s\n',e.message)
%end

fprintf('\nPP_Completed\n');

out = 1;
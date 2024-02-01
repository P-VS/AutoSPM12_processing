function out = my_spmbatch_run_aslpreprocessing(sub,ses,datpath,paramsfile)

load(paramsfile)

try
    %% preprocess anatomical scans
    if params.asl.do_segmentation
        % Normalization
        params.vbm.do_normalization = true;
        params.vbm.normvox = params.asl.normvox;
    
        % Segmentation
        params.vbm.do_segmentation = true;
        params.vbm.do_roi_atlas = false;
        params.vbm.do_surface = false;
    
        [delfiles,keepfiles] = my_spmbatch_cat12vbm(sub,ses,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params.save_intermediate_results,'anat','preproc_anat_asl');
    end
    
    %% preprocess asl scans
    if params.asl.preprocessed_ge
        [delfiles,keepfiles] = my_spmbatch_aslpreprocessed(sub,ses,datpath,params);
    
        % Clean up unnecessary files
        cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,params.save_intermediate_results,'perf',params.save_folder);
    end
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

fprintf('\nPP_Completed\n');

out = 1;
function out = my_spmbatch_run_fmriprocessing(sub,ses,run,task,datpath,paramsfile)

load(paramsfile)

try
    %% make batch
    matlabbatch = my_spmbatch_fmrilevel1processing(sub,ses,run,task,datpath,params);

    if ~isempty(matlabbatch), spm_jobman('run', matlabbatch); end
catch e
    fprintf('\nPP_Error\n');
    fprintf('\nThe error was: \n%s\n',e.message)
end

fprintf('\nPP_Completed\n');

out = 1;
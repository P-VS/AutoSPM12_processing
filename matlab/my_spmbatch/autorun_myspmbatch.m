function [rpfile,ppfile] = autorun_myspmbatch(sub,ses,task,datpath,fieldmap,pepolar,nechoes,dummytime,do_realignment,do_slicetime,do_segmentation,do_normalization,do_smoothing)

[matlabbatch,delfiles,keepfiles] = my_spmbatch(sub,ses,task,datpath,nechoes,dummytime,fieldmap,pepolar,do_realignment,do_slicetime,do_segmentation,do_normalization,do_smoothing);

%% Run batch
fprintf(['\nStart preprocessing'])

spm_jobman('run', matlabbatch);

cleanup_intermediate_files(sublist(i),nsessions(j),datpath,delfiles,keepfiles,save_intermediate_results)
   
fprintf(['\nDone\n'])

if sub<10
    substring = ['sub-0' num2str(sub)];
else
    substring = ['sub-' num2str(sub)];
end

subpath = fullfile(datpath,substring,c);
preproc_func = fullfile(subpath,'preproc_func');

rpdir = dir([preproc_func filesep 'rp_*' task '_bold.txt']);
rpfile = fullfile(preproc_func,rpdir.name);

ppdir = dir([preproc_func filesep 'w*' task '_bold.nii']);
ppfile = fullfile(preproc_func,ppdir.name);

end


%%_________________________________________________________________________________________
function cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,save_intermediate_results)

if sub<10
    substring = ['sub-0' num2str(sub)];
else
    substring = ['sub-' num2str(sub)];
end

subpath = fullfile(datpath,substring,['ses-00' num2str(ses)]);

subanatdir = fullfile(subpath,'anat');
subfmridir = fullfile(subpath,'func');

preproc_anat = fullfile(subpath,'preproc_anat_test');
preproc_func = fullfile(subpath,'preproc_func_test');

if ~exist(preproc_anat, 'dir')
    mkdir(preproc_anat);
end
if ~exist(preproc_func, 'dir')
    mkdir(preproc_func);
end


%% Delete intermediate files
if ~save_intermediate_results
    for i=1:numel(delfiles)
        delete(delfiles{i}{1});
    end
end

%% Move results
for i=1:numel(keepfiles)
    if contains(keepfiles{i},'anat')
        movefile(keepfiles{i}{1},preproc_anat)
    end
    if contains(keepfiles{i},'func')
        movefile(keepfiles{i}{1},preproc_func)
    end
end

end
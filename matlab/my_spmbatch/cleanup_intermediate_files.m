function cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,save_intermediate_results,params)

if sub<10
    substring = ['sub-0' num2str(sub)];
else
    substring = ['sub-' num2str(sub)];
end

subpath = fullfile(datpath,substring,['ses-00' num2str(ses)]);

subanatdir = fullfile(subpath,'anat');
subfmridir = fullfile(subpath,'func');

preproc_anat = fullfile(subpath,'preproc_anat');
preproc_func = fullfile(subpath,params.save_folder);

if ~exist(preproc_anat, 'dir')
    mkdir(preproc_anat);
end
if ~exist(preproc_func, 'dir')
    mkdir(preproc_func);
end


%% Delete intermediate files
if ~save_intermediate_results
    for i=1:numel(delfiles)
        if isfile(delfiles{i})
            if contains(delfiles{i},'anat')
                if isfile(delfiles{i}{1}); delete(delfiles{i}{1}); end
            end
            if contains(delfiles{i},'func')
                if isfile(delfiles{i}{1}); delete(delfiles{i}{1}); end
            end
            if contains(delfiles{i},'fmap')
                if isfile(delfiles{i}{1}); delete(delfiles{i}{1}); end
            end
        elseif isfolder(delfiles{i})
            rmdir(delfiles{i}{1},'s');
        end
    end
end

%% Move results
for i=1:numel(keepfiles)
    if contains(keepfiles{i},'anat')
        if isfile(keepfiles{i}{1})
            [opath,~,~] = fileparts(keepfiles{i}{1});
            if ~strcmp(opath,preproc_anat)
                movefile(keepfiles{i}{1},preproc_anat); 
            end
        end
    end
    if contains(keepfiles{i},'func')
        if isfile(keepfiles{i}{1}) 
            [opath,~,~] = fileparts(keepfiles{i}{1});
            if ~strcmp(opath,preproc_func)
                movefile(keepfiles{i}{1},preproc_func); 
            end
        end
    end
end

end
function cleanup_intermediate_files(sub,ses,datpath,delfiles,keepfiles,save_intermediate_results,subdir,save_folder)

if sub<10
    substring = ['sub-0' num2str(sub)];
else
    substring = ['sub-' num2str(sub)];
end

subpath = fullfile(datpath,substring,['ses-00' num2str(ses)]);

subolddir = fullfile(subpath,subdir);
subnewdir = fullfile(subpath,save_folder);

if ~exist(subnewdir, 'dir')
    mkdir(subnewdir);
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
    if isfile(keepfiles{i}{1})
        [opath,~,~] = fileparts(keepfiles{i}{1});
        if ~strcmp(opath,subnewdir)
            movefile(keepfiles{i}{1},subnewdir); 
        end
    end
end

end
function my_spmbatch_asllevel1processing

matlabbatch = {};

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub,['%0' num2str(params.sub_digits) 'd'])];

ppparams.sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring,ppparams.sesstring)), ppparams.sesstring = ['ses-' num2str(ses,'%03d')]; end

ppparams.subpath = fullfile(datpath,ppparams.substring,ppparams.sesstring);

if ~isfolder(ppparams.subpath), ppparams.subpath = fullfile(datpath,ppparams.substring); end

if ~isfolder(ppparams.subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.subfuncdir = fullfile(ppparams.subpath,'func');

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

ppparams.preprocfmridir = fullfile(ppparams.subpath,params.preprocfmridir);

if ~isfolder(ppparams.subfuncdir)
    fprintf(['No preprocessed func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

for ir=1:numel(params.iruns)
    %% Search for the data files
    
    namefilters(1).name = ppparams.substring;
    namefilters(1).required = true;
    
    namefilters(2).name = ppparams.sesstring;
    namefilters(2).required = false;
    
    switch params.func.use_runs
        case 'separately'
            namefilters(3).name = ['run-' num2str(run)];
        case 'together'
            namefilters(3).name = ['run-' num2str(params.iruns(ir))];
    end
    if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end
    
    namefilters(4).name = ['task-' task];
    namefilters(4).required = true;
    
    %fASL data
    
    if contains(params.modality,'fasl'), namefilters(5).name = '_cbf'; end
    namefilters(5).required = true;

    namefilters(6).name = params.fmri_prefix;
    namefilters(6).required = true;

    funcniilist = my_spmbatch_dirfilelist(ppparams.preprocfmridir,'nii',namefilters,false);
    
    if isempty(funcniilist)
        fprintf(['No nii files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' params.task '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    if params.func.meepi && params.use_echoes_as_sessions %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(params.func.echoes(1))]),contains({funcniilist.name},['_e' num2str(params.func.echoes(1))])));
        if isempty(tmp), edirniilist = funcniilist; else edirniilist = funcniilist(tmp); end
    else
        edirniilist = funcniilist;
    end

    prefixlist = split({edirniilist.name},'sub-');
    if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

    tmp = find(strcmp(prefixlist,params.fmri_prefix));
    if ~isempty(tmp), ppparams.frun(ir).func(1).funcfile = edirniilist(tmp).name; end

    if ~isfield(ppparams.frun(ir).func(1),'funcfile')
        fprintf(['no preprocessed fmri data found for run ' num2str(ir) ' for echo ' num2str(params.func.echoes(1)) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    Vfunc = spm_vol(fullfile(ppparams.preprocfmridir,ppparams.frun(ir).func(1).funcfile));

    for i=1:numel(Vfunc)
        ppfmridat{ir}.sess{1}.func{i,1} = [Vfunc(i).fname ',' num2str(i)];
    end
    
    %events.tsv file
    
    enamefilters(1:4) = namefilters(1:4);
    
    enamefilters(5).name = '_events';
    enamefilters(5).required = true;
    
    functsvlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'tsv',enamefilters,false);
    
    if isempty(functsvlist)
        fprintf(['No events.tsv files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
        fprintf('\nPP_Error\n');
        return
    end
    
    ppparams.frun(ir).functsvfile = fullfile(functsvlist(1).folder,functsvlist(1).name);
    
    %json file

    jnamefilters(1:4) = namefilters(1:4);
    
    jnamefilters(5).name = '_asl';
    jnamefilters(5).required = true;
    
    funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',jnamefilters,false);
        
    jstmp = find(or(contains({funcjsonlist.name},'echo-1'),contains({funcjsonlist.name},'_e1')));
    funcjsonlist = funcjsonlist(jstmp);
    
    if isempty(funcjsonlist)
        fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' task '\n'])
        fprintf('\nPP_Error\n');
        return
    end
    
    ppparams.frun(ir).funcjsonfile = fullfile(funcjsonlist(1).folder,funcjsonlist(1).name);
end

%% fMRI model specification
if params.func.mruns && contains(params.func.use_runs,'separately')
    resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname '_run-' num2str(run)]);
else
    resultmap = fullfile(ppparams.subpath,['SPMMAT-' task '_' params.analysisname]);
end

if exist(resultmap,'dir'); rmdir(resultmap,'s'); end
mkdir(resultmap)

%% Model Specification

jsondat = fileread(ppparams.frun(1).funcjsonfile);
jsondat = jsondecode(jsondat);

if ~params.reduced_temporal_resolution, tr = jsondat.RepetitionTime; else tr=params.newTR; end

for ir=1:numel(params.iruns)
    % correct events file for dummy scans if needed
    dummys = floor(params.dummytime/tr);
    try
        edat{ir} = tdfread(ppparams.frun(ir).functsvfile,'\t');
    catch
        T = readtable(ppparams.frun(ir).functsvfile,'FileType','text');
        edat{ir}.onset = T.onset;
        edat{ir}.duration = T.duration;
        edat{ir}.trial_type = T.trial_type;
    end
    edat{ir}.onset = edat{ir}.onset-dummys*tr;
    
    for it=1:numel(edat{ir}.trial_type(:,1))
        ntrial_type(it,1) = convertCharsToStrings(edat{ir}.trial_type(it,:));
    end
    
    edat{ir}.trial_type = ntrial_type;
    
    [~,edatorder] = sort(edat{ir}.trial_type);
    edat{ir}.onset = edat{ir}.onset(edatorder);
    edat{ir}.duration = edat{ir}.duration(edatorder);
    edat{ir}.trial_type = edat{ir}.trial_type(edatorder,:);

    numc=0;
    for trial=1:numel(edat{ir}.onset)
        if isstring(edat{ir}.trial_type(trial,:))
            trial_type = edat{ir}.trial_type(trial,:);
        else
            trial_type = num2str(edat{ir}.trial_type(trial,:));
        end
    
        if numc>0
            numc = numel(edat{ir}.conditions)+1;
            for nc=1:numel(edat{ir}.conditions)
                if strcmp(edat{ir}.conditions{nc}.name,strtrim(trial_type)); numc=nc; end
            end
            if numc<numel(edat{ir}.conditions)+1
                edat{ir}.conditions{numc}.onsets = [edat{ir}.conditions{numc}.onsets edat{ir}.onset(trial)];
                edat{ir}.conditions{numc}.durations = [edat{ir}.conditions{numc}.durations edat{ir}.duration(trial)];
            else
                edat{ir}.conditions{numc}.name = strtrim(trial_type);
                edat{ir}.conditions{numc}.onsets = [edat{ir}.onset(trial)];
                edat{ir}.conditions{numc}.durations = [edat{ir}.duration(trial)];
            end
        else
            edat{ir}.conditions{1}.name = strtrim(trial_type);
            edat{ir}.conditions{1}.onsets = edat{ir}.onset(trial);
            edat{ir}.conditions{1}.durations = edat{ir}.duration(trial);
            numc=1;
        end
    end
    
    for nc=1:numel(edat{1}.conditions)
        edat{ir}.conditions{nc}.startblock=floor(1+edat{ir}.conditions{nc}.onsets/tr);
        edat{ir}.conditions{nc}.endblock=floor((edat{ir}.conditions{nc}.onsets+edat{ir}.conditions{nc}.durations)/tr);
    end
    
end
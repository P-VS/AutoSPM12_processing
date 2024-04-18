function [params] = my_spmbatch_fmriprocessing_do_smoothing(sub,ses,run,task,datpath,params)

%% Search for the data folders

substring = ['sub-' num2str(sub,'%02d')];
if ~isfolder(fullfile(datpath,substring)), substring = ['sub-' num2str(sub,'%03d')]; end

sesstring = ['ses-' num2str(ses,'%02d')];
if ~isfolder(fullfile(datpath,substring,sesstring)), sesstring = ['ses-' num2str(ses,'%03d')]; end

subpath = fullfile(datpath,substring,sesstring);

if ~isfolder(subpath), subpath = fullfile(datpath,substring); end

if ~isfolder(subpath)
    fprintf(['No data folder for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

subfuncdir = fullfile(subpath,params.preprocfmridir);

if ~isfolder(subfuncdir)
    fprintf(['No func data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
end

%% Search for the data files

namefilters(1).name = substring;
namefilters(1).required = true;

namefilters(2).name = sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(run)];
if params.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['task-' task];
namefilters(4).required = true;

namefilters(5).name = '_bold';
namefilters(5).required = true;

namefilters(6).name = params.fmri_prefix;
namefilters(6).required = true;

funcniilist = my_spmbatch_dirfilelist(subfuncdir,'nii',namefilters,false);

if isempty(funcniilist)
    fprintf(['No nifti files found for ' substring ' ' sesstring ' task-' task '\n'])
    fprintf('\nPP_Error\n');
    return
end

%% do smoothing
for ie=params.echoes
    if params.meepi %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(ie)]),contains({funcniilist.name},['_e' num2str(ie)])));
        if isempty(tmp), funcfile = funcniilist(1).name; else funcfile = funcniilist(tmp).name; end
    else
        funcfile = funcniilist(1).name;
    end

    fprintf('Do smoothing \n')

    Vfunc = spm_vol(fullfile(subfuncdir,funcfile));

    wfuncdat = spm_read_vols(Vfunc);

    Vout = Vfunc;
    swfuncdat = zeros(size(wfuncdat));

    spm_progress_bar('Init',numel(Vfunc),'Smoothing','volumes completed');

    for i=1:numel(Vfunc)
        tswfuncdat = my_spmbatch_smooth(wfuncdat(:,:,:,i),Vfunc(i),[],[6 6 6],0);

        swfuncdat(:,:,:,i) = tswfuncdat;

        spm_progress_bar('Set',i);
    end

    spm_progress_bar('Clear');

    [pth,nm,~] = fileparts(Vfunc(i).fname);
    for j=1:numel(Vout)
        Vout(j).fname = fullfile(pth, ['s' nm  '.nii']);
        Vout(j).descrip = 'my_spmbatch - smooth';
        Vout(j).pinfo = [1,0,0];
        Vout(j).n = [j 1];
    end

    Vout = myspm_write_vol_4d(Vout,swfuncdat);

    fprintf('Done smoothing \n')

    clear wfuncdat
end

params.fmri_prefix = ['s' params.fmri_prefix];
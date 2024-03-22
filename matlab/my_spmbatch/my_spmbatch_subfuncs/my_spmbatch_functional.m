function [delfiles,keepfiles] = my_spmbatch_functional(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

if ~params.func.meepi, params.func.echoes = [1]; end

ppparams.sub = sub;
ppparams.ses = ses;
ppparams.run = run;
ppparams.task = task;
ppparams.reorient = params.reorient;
ppparams.use_parallel = params.use_parallel;

%% Search for the data folders

ppparams.substring = ['sub-' num2str(sub,'%02d')];
if ~isfolder(fullfile(datpath,ppparams.substring)), ppparams.substring = ['sub-' num2str(sub,'%03d')]; end

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

if params.func.fieldmap || params.func.pepolar
    ppparams.subfmapdir = fullfile(ppparams.subpath,'fmap');

    if ~isfolder(ppparams.subfmapdir)
        fprintf(['No fmap data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
        fprintf('\nPP_Error\n');

        ppparams.subfmapdir = ppparams.subfuncdir;
    end
end

%% Search for the data files

namefilters(1).name = ppparams.substring;
namefilters(1).required = true;

namefilters(2).name = ppparams.sesstring;
namefilters(2).required = false;

namefilters(3).name = ['run-' num2str(ppparams.run)];
if params.func.mruns, namefilters(3).required = true; else namefilters(3).required = false; end

namefilters(4).name = ['task-' ppparams.task];
namefilters(4).required = true;

namefilters(5).name = '_bold';
namefilters(5).required = true;

funcniilist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'nii',namefilters,false);
funcjsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',namefilters,false);

if isempty(funcniilist)
    fprintf(['No nifti files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

if isempty(funcjsonlist)
    fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

if params.func.fieldmap || params.func.pepolar
    
    fmfilters(1).name = ppparams.substring;
    fmfilters(1).required = true;
    
    fmfilters(2).name = ppparams.sesstring;
    fmfilters(2).required = false;
    
    fmfilters(3).name = ['run-' num2str(ppparams.run)];
    fmfilters(3).required = false;

    if params.func.pepolar  
        fmfilters(4).name = 'dir-';
        fmfilters(4).required = true;

        fmfilters(4).name = '_epi';
        fmfilters(4).required = true;
    
        fmapniilist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'nii',fmfilters,true);
        fmapjsonlist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'json',fmfilters,true);
    
        if isempty(fmapniilist) || isempty(fmapjsonlist)
            params.func.pepolar = false;
    
            fprintf(['No fmap files found for ' ppparams.substring ' ' ppparams.sesstring ' run ' num2str(ppparams.run) ' task-' ppparams.task '\n'])
            fprintf('\nPP_Error\n');
        end
    elseif params.func.fieldmap

        fmfilters(4).name = '_phase';
        fmfilters(4).required = true;
    
        phniilist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'nii',fmfilters,true);
    
        fmfilters(4).name = '_magnitude';
        fmfilters(4).required = true;
    
        amniilist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'nii',fmfilters,true);
        fmapjsonlist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'json',fmfilters,true);

        if isempty(phniilist) || isempty(amniilist) || isempty(fmapjsonlist)
            params.func.fieldmap = false;
        else
            tmp = find(contains({phniilist.name},'_phase1'));
            if isempty(tmp), params.func.fieldmap = false; else ppparams.fmap.phim1 = fullfile(phniilist(tmp).folder,phniilist(tmp).name); end

            tmp = find(contains({phniilist.name},'_phase2'));
            if isempty(tmp), params.func.fieldmap = false; else ppparams.fmap.phim2 = fullfile(phniilist(tmp).folder,phniilist(tmp).name); end

            tmp = find(contains({amniilist.name},'_magnitude1'));
            if isempty(tmp), params.func.fieldmap = false; else ppparams.fmap.amim1 = fullfile(amniilist(tmp).folder,amniilist(tmp).name); end

            tmp = find(contains({amniilist.name},'_magnitude2'));
            if isempty(tmp), params.func.fieldmap = false; else ppparams.fmap.amim2 = fullfile(amniilist(tmp).folder,amniilist(tmp).name); end

            tmp = find(contains({fmapjsonlist.name},'_magnitude1'));
            if isempty(tmp), params.func.fieldmap = false; else ppparams.fmap.json1 = fullfile(fmapjsonlist(tmp).folder,fmapjsonlist(tmp).name); end

            tmp = find(contains({fmapjsonlist.name},'_magnitude2'));
            if isempty(tmp), params.func.fieldmap = false; else ppparams.fmap.json2 = fullfile(fmapjsonlist(tmp).folder,fmapjsonlist(tmp).name); end
        end
    end
end

for ie=params.func.echoes
    if params.func.meepi %Filter list based on echo number
        tmp = find(or(contains({funcniilist.name},['_echo-' num2str(ie)]),contains({funcniilist.name},['_e' num2str(ie)])));
        if isempty(tmp)
            fprintf(['no fmri data found for echo ' num2str(ie) '\n'])
            fprintf('\nPP_Error\n');
            return
        end
    
        edirniilist = funcniilist(tmp);

        if ie==params.func.echoes(1)
            tmp = find(and(~contains({funcniilist.name},'_echo-'),~contains({funcniilist.name},'_e')));
            if ~isempty(tmp)
                edirniilist = [edirniilist;funcniilist(tmp)];
            end
        end

        jstmp = find(or(contains({funcjsonlist.name},['echo-' num2str(ie)]),contains({funcjsonlist.name},['_e' num2str(ie)])));
        if isempty(jstmp)
            fprintf(['no json file found for echo ' num2str(ie) '\n'])
            fprintf('\nPP_Error\n');
            return
        end

        ppparams.func(ie).jsonfile = fullfile(funcjsonlist(jstmp(1)).folder,funcjsonlist(jstmp(1)).name);
    else
        edirniilist = funcniilist;

        ppparams.func(ie).jsonfile = fullfile(funcjsonlist(1).folder,funcjsonlist(1).name);
    end

    prefixlist = split({edirniilist.name},'sub-');
    if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

    tmp = find(strlength(prefixlist)==0);
    if ~isempty(tmp), ppparams.func(ie).funcfile = edirniilist(tmp).name; end

    ppparams.func(ie).prefix = '';
    studyprefix = 'e';

    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.func(ie).efuncfile = edirniilist(tmp).name;
        ppparams.func(ie).prefix = studyprefix;

        tmp = find(strcmp(prefixlist,'fe'));
        if ~isempty(tmp),
            ppparams.func(ie).fefuncfile = edirniilist(tmp).name;
            ppparams.func(ie).reffunc = ppparams.func(ie).fefuncfile;
        end
    end

    if params.func.do_realignment, studyprefix = ['r' studyprefix]; end

    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.func(ie).rfuncfile = edirniilist(tmp).name;

        tmp = find(strcmp(prefixlist,['mean' ppparams.func(ie).prefix]));
        if ~isempty(tmp)
            ppparams.func(ie).meanfuncfile = edirniilist(tmp(1)).name;
            ppparams.func(ie).reffunc = ppparams.func(ie).meanfuncfile;
        end

        ppparams.func(ie).prefix = studyprefix;
    end

    if params.func.fieldmap || params.func.pepolar, studyprefix = ['u' studyprefix]; end

    tmp = find(strcmp(prefixlist,[ppparams.func(ie).prefix studyprefix]));
    if ~isempty(tmp)
        ppparams.func(ie).ufuncfile = edirniilist(tmp).name;

        tmp = find(strcmp(prefixlist,['fu' ppparams.func(ie).prefix]));
        if ~isempty(tmp), ppparams.func(ie).fefuncfile = edirniilist(tmp).name; end

        ppparams.func(ie).reffunc = ppparams.func(ie).fefuncfile;

        tmp = find(strcmp(prefixlist,['umean' ppparams.func(ie).prefix]));
        if ~isempty(tmp), ppparams.func(ie).meanfuncfile = edirniilist(tmp).name; end

        ppparams.func(ie).reffunc = ppparams.func(ie).meanfuncfile;

        ppparams.func(ie).prefix = studyprefix;
    end

    if params.func.do_slicetime, studyprefix = ['a' studyprefix]; end

    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.func(ie).afuncfile = edirniilist(tmp).name;
        ppparams.func(ie).prefix = studyprefix;
    end

    if params.func.do_echocombination, studyprefix = ['c' studyprefix]; end

    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.func(ie).cfuncfile = edirniilist(tmp).name;
        ppparams.func(ie).prefix = studyprefix;
    end

    if params.func.do_normalization, studyprefix = ['w' studyprefix]; end
    
    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.func(ie).wfuncfile = edirniilist(tmp).name;
        ppparams.func(ie).prefix = studyprefix;
    end

    if params.func.do_smoothing, studyprefix = ['s' studyprefix]; end

    tmp = find(strcmp(prefixlist,studyprefix));
    if ~isempty(tmp)
        ppparams.func(ie).sfuncfile = edirniilist(tmp).name;
        ppparams.func(ie).prefix = studyprefix;
    end

    if params.func.pepolar   
        if params.func.meepi %Filter list based on echo number
            tmp = find(or(contains({fmapniilist.name},['echo-' num2str(ie)]),contains({fmapniilist.name},['_e' num2str(ie)])));
            if isempty(tmp)
                fprintf(['no fmap data found for echo ' num2str(ie) '\n'])
                return
            end
        
            efmapniilist = fmapniilist(tmp);
    
            jstmp = find(or(contains({fmapjsonlist.name},['echo-' num2str(ie)]),contains({fmapjsonlist.name},['_e' num2str(ie)])));
            if isempty(jstmp)
                fprintf(['no fmap json file found for echo ' num2str(ie) '\n'])
                fprintf('\nPP_Error\n');
                return
            end
    
            ppparams.func(ie).fmapjsonfile = fullfile(fmapjsonlist(jstmp(1)).folder,fmapjsonlist(jstmp(1)).name);

        else
            efmapniilist = fmapniilist;
    
            ppparams.func(ie).fmapjsonfile = fullfile(fmapjsonlist(1).folder,fmapjsonlist(1).name);
        end
    
        fmprefixlist = split({efmapniilist.name},'sub-');
        fmprefixlist = fmprefixlist(:,:,1);
            
        tmp = find(strlength(fmprefixlist)==0);
        if ~isempty(tmp)
            ppparams.func(ie).fmapfile = efmapniilist(tmp).name; 
        else
            fprintf(['no fmap data found for echo ' num2str(ie) '\n'])
            fprintf('\nPP_Error\n');
            params.func.pepolar = false;
        end

        ppparams.func(ie).fmap_prefix = '';
    end
end

%% Do the preprocessing
[delfiles,keepfiles] = my_spmbatch_preprocfunc(ppparams,params,delfiles,keepfiles);
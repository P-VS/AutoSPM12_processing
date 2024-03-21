function [delfiles,keepfiles] = my_spmbatch_fasl(sub,ses,run,task,datpath,params)

delfiles = {};
keepfiles = {};

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

ppparams.subperfdir = fullfile(ppparams.subpath,'perf');

if ~isfolder(ppparams.subperfdir)
    fprintf(['No perf data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
    fprintf('\nPP_Error\n');
    return
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

% asl data

namefilters(5).name = '_asl';
namefilters(5).required = true;

aslniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if isempty(aslniilist)
    fprintf(['No asl nifti files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

% m0scan data

namefilters(5).name = '_m0scan';
namefilters(5).required = true;

m0scaniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if isempty(m0scaniilist)
    fprintf(['No m0scan nifti files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

% json data

namefilters(5).name = '_bold';
namefilters(5).required = true;

asljsonlist = my_spmbatch_dirfilelist(ppparams.subperfdir,'json',namefilters,false);

if isempty(asljsonlist), asljsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',namefilters,false); end

if isempty(asljsonlist)
    fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

%% Search for the data files per echo

for ie=params.func.echoes
    % asl data files
    tmp = find(or(contains({aslniilist.name},['_echo-' num2str(ie)]),contains({aslniilist.name},['_e' num2str(ie)])));
    if isempty(tmp)
        fprintf(['no asl data found for echo ' num2str(ie) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    easlniilist = aslniilist(tmp);

    aslprefixlist = split({easlniilist.name},'sub-');
    if numel(easlniilist)==1, aslprefixlist=aslprefixlist{1}; else aslprefixlist = aslprefixlist(:,:,1); end

    tmp = find(strlength(aslprefixlist)==0);
    if ~isempty(tmp), ppparams.asl(ie).aslfile = easlniilist(tmp).name; end

    ppparams.asl(ie).aslprefix = '';
    studyprefix = 'e';

    if ~isempty(find(contains(aslprefixlist,'e')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).easlfile = easlniilist(tmp).name;
        end
        ppparams.asl(ie).aslprefix = studyprefix;
    end

    studyprefix = ['r' ppparams.asl(ie).aslprefix];
    if ~isempty(find(contains(aslprefixlist,'r')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).raslfile = easlniilist(tmp).name;
        end
        ppparams.asl(ie).aslprefix = studyprefix;
    end

    studyprefix = ['u' ppparams.asl(ie).aslprefix];
    if ~isempty(find(contains(aslprefixlist,'u')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).uaslfile = easlniilist(tmp).name;
        end
        ppparams.asl(ie).aslprefix = studyprefix;
    end

    % m0scan data files
    tmp = find(or(contains({m0scaniilist.name},['echo-' num2str(ie)]),contains({m0scaniilist.name},['_e' num2str(ie)])));
    if isempty(tmp)
        fprintf(['no m0scan file found for echo ' num2str(ie) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    em0scaniilist = m0scaniilist(tmp);

    m0prefixlist = split({em0scaniilist.name},'sub-');
    if numel(em0scaniilist)==1, m0prefixlist=m0prefixlist{1}; else m0prefixlist = m0prefixlist(:,:,1); end

    tmp = find(strlength(m0prefixlist)==0);
    if ~isempty(tmp), ppparams.asl(ie).m0scanfile = em0scaniilist(tmp).name; end

    ppparams.asl(ie).m0scanprefix = '';
    studyprefix = 'e';

    if ~isempty(find(contains(m0prefixlist,'e')))
        tmp = find(strcmp(m0prefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).em0scanfile = em0scaniilist(tmp).name;
        end
        ppparams.asl(ie).m0scanprefix = studyprefix;
    end

    studyprefix = ['r' ppparams.asl(ie).m0scanprefix];
    if ~isempty(find(contains(m0prefixlist,'r')))
        tmp = find(strcmp(m0prefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).rm0scanfile = em0scaniilist(tmp).name;
        end
        ppparams.asl(ie).m0scanprefix = studyprefix;
    end

    studyprefix = ['u' ppparams.asl(ie).m0scanprefix];
    if ~isempty(find(contains(m0prefixlist,'u')))
        tmp = find(strcmp(m0prefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).um0scanfile = em0scaniilist(tmp).name;
        end
        ppparams.asl(ie).m0scanprefix = studyprefix;
    end

    % json files
    jstmp = find(or(contains({asljsonlist.name},['echo-' num2str(ie)]),contains({asljsonlist.name},['_e' num2str(ie)])));
    if isempty(jstmp)
        fprintf(['no json file found for echo ' num2str(ie) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    ppparams.asl(ie).asljsonlist = fullfile(asljsonlist(jstmp(1)).folder,asljsonlist(jstmp(1)).name);

    if (~isfield(ppparams.asl(ie),'rm0scanfile') && ~isfield(ppparams.asl(ie),'um0scanfile')) ...
            || (~isfield(ppparams.asl(ie),'raslfile') && ~isfield(ppparams.asl(ie),'uaslfile'))
    end
end

%% do the preprocessing

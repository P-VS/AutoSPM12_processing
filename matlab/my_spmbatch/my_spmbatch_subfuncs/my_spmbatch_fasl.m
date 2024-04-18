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

% json data

namefilters(5).name = '_asl';
namefilters(5).required = true;

asljsonlist = my_spmbatch_dirfilelist(ppparams.subperfdir,'json',namefilters,false);

if isempty(asljsonlist)
    namefilters(5).name = '_bold';
    namefilters(5).required = true;

    asljsonlist = my_spmbatch_dirfilelist(ppparams.subfuncdir,'json',namefilters,false); 
end

if isempty(asljsonlist)
    fprintf(['No json files found for ' ppparams.substring ' ' ppparams.sesstring ' task-' ppparams.task '\n'])
    fprintf('\nPP_Error\n');
    return
end

% rp_ file

namefilters(5).name = 'rp_';
namefilters(5).required = true;

rplist = my_spmbatch_dirfilelist(ppparams.subperfdir,'.txt',namefilters,false);

if ~isempty(rplist), ppparams.asl(1).rp_file = fullfile(rplist(1).folder,rplist(1).name); end

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

    tmp = find(strcmp(aslprefixlist,'fe'));
    if ~isempty(tmp)
        ppparams.asl(ie).feaslfile = easlniilist(tmp).name;
        ppparams.asl(ie).refasl = ppparams.asl(ie).feaslfile;

        tmp = find(~strcmp(aslprefixlist,'fe'));
        aslprefixlist = aslprefixlist(tmp);
        easlniilist = easlniilist(tmp);
    end

    tmp = find(strcmp(aslprefixlist,'meane'));
    if ~isempty(tmp)
        ppparams.asl(ie).meanaslfile = easlniilist(tmp).name;
        ppparams.asl(ie).refasl = ppparams.asl(ie).meanaslfile;

        tmp = find(~strcmp(aslprefixlist,'meane'));
        aslprefixlist = aslprefixlist(tmp);
        easlniilist = easlniilist(tmp);
    end

    ppparams.asl(ie).aslprefix = '';
    studyprefix = 'e';

    if ~isempty(find(contains(aslprefixlist,'e')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).easlfile = easlniilist(tmp).name;
    
            if ~isfield(ppparams.asl(ie),'aslfile') || isempty(ppparams.asl(ie).aslfile)
                splitaslfile = split(ppparams.asl(ie).easlfile,'sub-');
                ppparams.asl(ie).aslfile = ['sub-' splitaslfile{2}];
            end

            ppparams.asl(ie).aslprefix = studyprefix;
        end
    end

    if params.func.do_realignment; studyprefix = ['r' studyprefix]; end
    
    if params.func.do_realignment && ~isempty(find(contains(aslprefixlist,'r')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).raslfile = easlniilist(tmp).name;

            if ~isfield(ppparams.asl(ie),'aslfile') || isempty(ppparams.asl(ie).aslfile)
                splitaslfile = split(ppparams.asl(ie).raslfile,'sub-');
                ppparams.asl(ie).aslfile = ['sub-' splitaslfile{2}];
            end

            ppparams.asl(ie).aslprefix = studyprefix;
        end
    end

    if params.func.pepolar; studyprefix = ['u' studyprefix]; end

    if params.func.pepolar && ~isempty(find(contains(aslprefixlist,'u')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).uaslfile = easlniilist(tmp).name;

            if ~isfield(ppparams.asl(ie),'aslfile') || isempty(ppparams.asl(ie).aslfile)
                splitaslfile = split(ppparams.asl(ie).uaslfile,'sub-');
                ppparams.asl(ie).aslfile = ['sub-' splitaslfile{2}];
            end

            ppparams.asl(ie).aslprefix = studyprefix;
        end
    end

    if params.asl.do_denoise; studyprefix = ['d' studyprefix]; end

    if params.asl.do_denoise && ~isempty(find(contains(aslprefixlist,'d')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(ie).daslfile = easlniilist(tmp).name;

            if ~isfield(ppparams.asl(ie),'aslfile') || isempty(ppparams.asl(ie).aslfile)
                splitaslfile = split(ppparams.asl(ie).daslfile,'sub-');
                ppparams.asl(ie).aslfile = ['sub-' splitaslfile{2}];
            end

            ppparams.asl(ie).aslprefix = studyprefix;
        end
    end

    if ~isfield(ppparams.asl(ie),'aslfile')
        fprintf(['no asl data found for echo ' num2str(ie) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    % m0scan data files
    if ~isempty(m0scaniilist)
        tmp = find(or(contains({m0scaniilist.name},['echo-' num2str(ie)]),contains({m0scaniilist.name},['_e' num2str(ie)])));
        if ~isempty(tmp)
    
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

                    if ~isfield(ppparams.asl(ie),'m0scanfile') || isempty(ppparams.asl(ie).m0scanfile)
                        splitaslfile = split(ppparams.asl(ie).em0scanfile,'sub-');
                        ppparams.asl(ie).m0scanfile = ['sub-' splitaslfile{2}];
                    end

                    ppparams.asl(ie).m0scanprefix = studyprefix;
                end
            end
        
            studyprefix = ['r' studyprefix];
            if ~isempty(find(contains(m0prefixlist,'r')))
                tmp = find(strcmp(m0prefixlist,studyprefix));
                if ~isempty(tmp)
                    ppparams.asl(ie).rm0scanfile = em0scaniilist(tmp).name;

                    if ~isfield(ppparams.asl(ie),'m0scanfile') || isempty(ppparams.asl(ie).m0scanfile)
                        splitaslfile = split(ppparams.asl(ie).rm0scanfile,'sub-');
                        ppparams.asl(ie).m0scanfile = ['sub-' splitaslfile{2}];
                    end

                    ppparams.asl(ie).m0scanprefix = studyprefix;
                end
            end
        
            if params.func.pepolar, studyprefix = ['u' studyprefix]; end

            if params.func.pepolar && ~isempty(find(contains(m0prefixlist,'u')))
                tmp = find(strcmp(m0prefixlist,studyprefix));
                if ~isempty(tmp)
                    ppparams.asl(ie).um0scanfile = em0scaniilist(tmp).name;

                    if ~isfield(ppparams.asl(ie),'m0scanfile') || isempty(ppparams.asl(ie).m0scanfile)
                        splitaslfile = split(ppparams.asl(ie).um0scanfile,'sub-');
                        ppparams.asl(ie).m0scanfile = ['sub-' splitaslfile{2}];
                    end

                    ppparams.asl(ie).m0scanprefix = studyprefix;
                end
            end
        end
    end

    % json files
    jstmp = find(or(contains({asljsonlist.name},['echo-' num2str(ie)]),contains({asljsonlist.name},['_e' num2str(ie)])));
    if isempty(jstmp)
        fprintf(['no json file found for echo ' num2str(ie) '\n'])
        fprintf('\nPP_Error\n');
        return
    end

    ppparams.asl(ie).asljson = fullfile(asljsonlist(jstmp(1)).folder,asljsonlist(jstmp(1)).name);
end

%% Data files after BOLD removal
% c..ASL data files
if ~isempty(aslniilist)
    tmp = find(and(~contains({aslniilist.name},'_echo-'),~contains({aslniilist.name},'_e')));
    if ~isempty(tmp)
        easlniilist = aslniilist(tmp);
    end
    
    aslprefixlist = split({easlniilist.name},'sub-');
    if numel(easlniilist)==1, aslprefixlist=aslprefixlist{1}; else aslprefixlist = aslprefixlist(:,:,1); end
    
    studyprefix = 'c';
    if params.asl.do_denoise; studyprefix = [studyprefix 'd']; end
    if params.func.pepolar; studyprefix = [studyprefix 'u']; end
    studyprefix = [studyprefix ppparams.asl(1).aslprefix];
    
    if params.asl.do_removebold && ~isempty(find(contains(aslprefixlist,'c')))
        tmp = find(strcmp(aslprefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(1).caslfile = aslniilist(tmp).name;
    
            splitaslfile = split(ppparams.asl(1).caslfile,'sub-');
            ppparams.asl(1).aslfile = ['sub-' splitaslfile{2}];
    
            ppparams.asl(1).aslprefix = studyprefix;
    
            params.func.echoes = [1];
        end
    end
end

% c..m0scan data files
if ~isempty(m0scaniilist)
    tmp = find(and(~contains({m0scaniilist.name},'_echo-'),~contains({m0scaniilist.name},'_e')));
    if ~isempty(tmp)
        em0scaniilist = m0scaniilist(tmp);
    end
    
    m0prefixlist = split({em0scaniilist.name},'sub-');
    if numel(m0prefixlist)==1, m0prefixlist=m0prefixlist{1}; else m0prefixlist = m0prefixlist(:,:,1); end
    
    studyprefix = 'c';
    if params.func.pepolar; studyprefix = [studyprefix 'u']; end
    studyprefix = [studyprefix ppparams.asl(1).m0scanprefix];
    
    if params.asl.do_removebold && ~isempty(find(contains(m0prefixlist,'c')))
        tmp = find(strcmp(m0prefixlist,studyprefix));
        if ~isempty(tmp)
            ppparams.asl(1).cm0scanfile = em0scaniilist(tmp).name;
    
            splitaslfile = split(ppparams.asl(1).cm0scanfile,'sub-');
            ppparams.asl(1).m0scanfile = ['sub-' splitaslfile{2}];
    
            ppparams.asl(1).m0scanprefix = studyprefix;
        end
    end
    
    if ~isempty(find(contains(m0prefixlist,'c1'))) || ~isempty(find(contains(m0prefixlist,'c2'))) || ~isempty(find(contains(m0prefixlist,'c2')))
        tmp = find(strcmp(m0prefixlist,'c1'));
        if ~isempty(tmp), ppparams.asl(ie).c1m0scanfile = em0scaniilist(tmp).name; end
    
        tmp = find(strcmp(m0prefixlist,'c2'));
        if ~isempty(tmp), ppparams.asl(ie).c2m0scanfile = em0scaniilist(tmp).name; end
    
        tmp = find(strcmp(m0prefixlist,'c3'));
        if ~isempty(tmp), ppparams.asl(ie).c3m0scanfile = em0scaniilist(tmp).name; end
    end
end

%% Search for existing CBF maps

namefilters(5).name = '_cbf';
namefilters(5).required = true;

cbfniilist = my_spmbatch_dirfilelist(ppparams.subperfdir,'nii',namefilters,false);

if ~isempty(cbfniilist)
    prefixlist = split({cbfniilist.name},'sub-');
    if numel(edirniilist)==1, prefixlist=prefixlist{1}; else prefixlist = prefixlist(:,:,1); end

    tmp = find(strlength(prefixlist)==0);
    if ~isempty(tmp), ppparams.cbfmap = fullfile(cbfniilist(tmp).folder,cbfniilist(tmp).name); end
end

if ~isfield(ppparams,'cbfmap'), params.asl.do_cbfmapping = true; else params.asl.do_cbfmapping = false; end

%% If geometric correction is needed
if params.func.pepolar && (~isfield(ppparams.asl(ie),'um0scanfile') || ~isfield(ppparams.asl(ie),'uaslfile'))
    ppparams.subfmapdir = fullfile(ppparams.subpath,'fmap');

    if ~isfolder(ppparams.subfmapdir)
        fprintf(['No fmap data folder found for subject ' num2str(sub) ' session ' num2str(ses)])
        fprintf('\nPP_Error\n');
    end

    fmfilters(1).name = ppparams.substring;
    fmfilters(1).required = true;
    
    fmfilters(2).name = ppparams.sesstring;
    fmfilters(2).required = false;
    
    fmfilters(3).name = ['run-' num2str(ppparams.run)];
    fmfilters(3).required = false;

    fmfilters(4).name = 'dir-';
    fmfilters(4).required = true;

    fmfilters(4).name = '_epi';
    fmfilters(4).required = true;

    fmapniilist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'nii',fmfilters,true);
    fmapjsonlist = my_spmbatch_dirfilelist(ppparams.subfmapdir,'json',fmfilters,true);

    if isempty(fmapniilist) || isempty(fmapjsonlist)
        params.func.pepolar = false;

        fprintf(['No fmap files found for ' ppparams.substring ' ' ppparams.sesstring ' run ' num2str(ppparams.run) ' task-' ppparams.task '\n'])
    end

    for ie=params.func.echoes
        if params.func.meepi %Filter list based on echo number
            tmp = find(or(contains({fmapniilist.name},['echo-' num2str(ie)]),contains({fmapniilist.name},['_e' num2str(ie)])));
            if isempty(tmp)
                fprintf(['no fmap data found for echo ' num2str(ie) '\n'])
                fprintf('\nPP_Error\n');
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
    end

    ppparams.func(ie).fmap_prefix = '';
end

%% Do fALS preprocessing
[delfiles,keepfiles] = my_spmbatch_preprocfasl(ppparams,params,delfiles,keepfiles);
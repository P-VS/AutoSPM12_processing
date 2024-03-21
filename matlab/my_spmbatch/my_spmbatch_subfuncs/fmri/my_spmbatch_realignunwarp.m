function [funcdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ne,ppparams,params,keepfiles,delfiles)

%-Realign
%--------------------------------------------------------------------------

if ne==params.func.echoes(1)
    eoptions.quality = 0.9;
    eoptions.sep = 4;
    eoptions.fwhm = 5;
    eoptions.rtm = 0;
    eoptions.interp = 2;
    eoptions.wrap = [0 0 0];
    eoptions.PW = '';

    if params.func.isaslbold
        spm_realign_asl(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile),eoptions);
    else
        spm_realign(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile),eoptions);
    end

    ppparams.func(1).rp_file = spm_file(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile), 'prefix','rp_','ext','.txt');

    if params.func.meepi
        [rppath,rpname,~] = fileparts(ppparams.func(1).rp_file);
        nrpname = split(rpname,'_echo-');
        nrp_file = fullfile(rppath,[nrpname{1} '_bold.txt']);

        movefile(ppparams.func(1).rp_file,nrp_file);

        ppparams.func(1).rp_file = nrp_file;
    end

    keepfiles{numel(keepfiles)+1} = {ppparams.func(1).rp_file};
    
end

fname = split(ppparams.func(ne).efuncfile,'.nii');
realign_mat = fullfile(ppparams.subfuncdir,[fname{1} '.mat']);

delfiles{numel(delfiles)+1} = {realign_mat};

if params.func.meepi && ne>params.func.echoes(1)
    Vtemp1 = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(params.func.echoes(1)).efuncfile));
    Vtemp2 = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile));
    for ti=1:numel(Vtemp2)
        spm_get_space([Vtemp2(ti).fname ',' num2str(ti)],Vtemp1(ti).mat);
    end
end

if params.func.fieldmap
    %% Realign and unwarp the func series
        
    jsondat = fileread(fullfile(ppparams.subfuncdir,ppparams.func(ne).jsonfile));
    jsondat = jsondecode(jsondat);
    pedir = jsondat.PhaseEncodingDirection;
    
    if contains(pedir,'i')
        pedim = 1;
        WrapD = [1 0 0];
    elseif contains(pedir,'j')
        pedim = 2;
        WrapD = [0 1 0];
    else
        pedim = 3;
        WrapD = [0 0 1];
    end
    
    %-Unwarp Estimate
    %--------------------------------------------------------------------------

    uweoptions.basfcn = [12 12];
    uweoptions.regorder = 1;
    uweoptions.lambda = 100000;
    uweoptions.jm = 0;
    uweoptions.fot = [4 5];
    uweoptions.sot = [];
    uweoptions.uwfwhm = 4;
    uweoptions.rem = 1;
    uweoptions.noi = 5;
    uweoptions.expround = 'Average';
    
    uweflags.sfP = ppparams.func(ne).vdm_file;
    P1 = deblank(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile));
    if isempty(spm_file(P1,'number'))
        P1 = spm_file(P1,'number',1);
    end
    VP1 = spm_vol(P1);
    uweflags.M = VP1.mat;
    ds = spm_uw_estimate(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile),uweoptions);
    sess(1).ds = ds;
    dsfile = spm_file(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile), 'suffix','_uw', 'ext','.mat');
    save(dsfile,'ds', spm_get_defaults('mat.format'));

    delfiles{numel(delfiles)+1} = {dsfile};
    
    %-Unwarp Write - Sessions should be within subjects
    %--------------------------------------------------------------------------

    uwroptions.uwwhich = [2 1];
    uwroptions.rinterp = 4;
    uwroptions.wrap = WrapD;
    uwroptions.mask = 1;
    uwroptions.prefix = 'u';
    
    funcdat = my_spmbatch_uw_apply(cat(2,sess.ds),uwroptions);

    ppparams.func(ne).meanfuncfile = ['meanu' ppparams.func(ne).efuncfile];
    ppparams.func(ne).reffunc = ppparams.func(ne).meanfuncfile;
    
    ppparams.func(ne).ufuncfile = ['u' ppparams.func(ne).efuncfile];
    ppparams.func(ne).prefix = ['u' ppparams.func(ne).prefix];

    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ne).meanfuncfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ne).ufuncfile)};
else
    %% Reslice the func series
    
    roptions.which = [2 1];
    roptions.interp = 4;
    roptions.wrap = [0 0 0];
    roptions.mask = 1;
    roptions.prefix = 'r';
    
    my_spmbatch_reslice(fullfile(ppparams.subfuncdir,ppparams.func(ne).efuncfile),roptions);

    ppparams.func(ne).meanfuncfile = ['mean' ppparams.func(ne).efuncfile];
    ppparams.func(ne).reffunc = ppparams.func(ne).meanfuncfile;

    ppparams.func(ne).rfuncfile = ['r' ppparams.func(ne).efuncfile];
    ppparams.func(ne).prefix = ['r' ppparams.func(ne).prefix];

    V = spm_vol(fullfile(ppparams.subfuncdir,ppparams.func(ne).rfuncfile));
    funcdat = spm_read_vols(V);

    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ne).meanfuncfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,ppparams.func(ne).rfuncfile)};
end
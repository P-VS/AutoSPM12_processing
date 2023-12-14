function [funcdat,ppparams,keepfiles,delfiles] = my_spmbatch_realignunwarp(ne,ppparams,params,keepfiles,delfiles)

%-Realign
%--------------------------------------------------------------------------

if ne==ppparams.echoes(1)
    eoptions.quality = 0.9;
    eoptions.sep = 4;
    eoptions.fwhm = 5;
    eoptions.rtm = 0;
    eoptions.interp = 2;
    eoptions.wrap = [0 0 0];
    eoptions.PW = '';

    spm_realign(ppparams.funcfile{ne},eoptions);

    ppparams.rp_file = spm_file(ppparams.funcfile{ne}, 'prefix','rp_','ext','.txt');

    if ppparams.meepi
        [rppath,rpname,~] = fileparts(ppparams.rp_file);
        nrpname = split(rpname,'bold_e');
        nrp_file = fullfile(rppath,[nrpname{1} 'bold.txt']);

        movefile(ppparams.rp_file,nrp_file);

        ppparams.rp_file = nrp_file;
    end

    keepfiles{numel(keepfiles)+1} = {ppparams.rp_file};
    
end

[fpath,fname,~] = fileparts(ppparams.funcfile{ne});
realign_mat = fullfile(fpath,[fname '.mat']);

delfiles{numel(delfiles)+1} = {realign_mat};

if ppparams.meepi && ne>ppparams.echoes(1)
    nrmname = split(ppparams.funcfile{ne},'bold_e');
    orealign_nii = [nrmname{1} 'bold_e' num2str(ppparams.echoes(1)) '.nii'];

    Vtemp1 = spm_vol(orealign_nii);
    Vtemp2 = spm_vol(ppparams.funcfile{ne});
    for ti=1:numel(Vtemp2)
        spm_get_space([Vtemp2(ti).fname ',' num2str(ti)],Vtemp1(ti).mat);
    end
end

if params.func.fieldmap
    %% Realign and unwarp the func series
        
    jsondat = fileread(ppparams.funcjsonfile);
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
    
    uweflags.sfP = ppparams.vdm_file;
    P1 = deblank(ppparams.funcfile{ne});
    if isempty(spm_file(P1,'number'))
        P1 = spm_file(P1,'number',1);
    end
    VP1 = spm_vol(P1);
    uweflags.M = VP1.mat;
    ds = spm_uw_estimate(ppparams.funcfile{ne},uweoptions);
    sess(1).ds = ds;
    dsfile = spm_file(ppparams.funcfile{ne}, 'suffix','_uw', 'ext','.mat');
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

    ppparams.reffunc{ne} = spm_file(ppparams.funcfile{ne}, 'prefix','meanu');
    ppparams.funcfile{ne} = spm_file(ppparams.funcfile{ne}, 'prefix','u');

    delfiles{numel(delfiles)+1} = {ppparams.reffunc{ne}};
    delfiles{numel(delfiles)+1} = {ppparams.funcfile{ne}};
else
    %% Reslice the func series
    
    roptions.which = [2 1];
    roptions.interp = 4;
    roptions.wrap = [0 0 0];
    roptions.mask = 1;
    roptions.prefix = 'r';
    
    my_spmbatch_reslice(ppparams.funcfile{ne},roptions);

    ppparams.reffunc{ne} = spm_file(ppparams.funcfile{ne}, 'prefix','mean');
    ppparams.funcfile{ne} = spm_file(ppparams.funcfile{ne}, 'prefix','r');

    V = spm_vol(ppparams.funcfile{ne});
    funcdat = spm_read_vols(V);

    delfiles{numel(delfiles)+1} = {ppparams.reffunc{ne}};
    delfiles{numel(delfiles)+1} = {ppparams.funcfile{ne}};
end
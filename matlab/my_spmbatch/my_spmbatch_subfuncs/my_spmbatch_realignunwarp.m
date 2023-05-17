function [reffunc,funcfile,rp_file,keepfiles,delfiles] = my_spmbatch_realignunwarp(subpath,substring,task,funcfile,vdm_file,params,keepfiles,delfiles)

if params.pepolar || params.fieldmap
    %% Realign and unwarp the func series

    if params.nechoes==1
        funcjsonfile = fullfile(subpath,'func',[substring '_task-' task '_bold.json']);
    else
        funcjsonfile = fullfile(subpath,'func',[substring '_task-' task '_bold_e1.json']);
    end
        
    jsondat = fileread(funcjsonfile);
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
    
    realignunwarp.data.scans = {funcfile};
    realignunwarp.data.pmscan(1) = {vdm_file};
    realignunwarp.eoptions.quality = 0.9;
    realignunwarp.eoptions.sep = 4;
    realignunwarp.eoptions.fwhm = 5;
    realignunwarp.eoptions.rtm = 0;
    realignunwarp.eoptions.einterp = 2;
    realignunwarp.eoptions.ewrap = [0 0 0];
    realignunwarp.eoptions.weight = '';
    realignunwarp.uweoptions.basfcn = [12 12];
    realignunwarp.uweoptions.regorder = 1;
    realignunwarp.uweoptions.lambda = 100000;
    realignunwarp.uweoptions.jm = 0;
    realignunwarp.uweoptions.fot = [4 5];
    realignunwarp.uweoptions.sot = [];
    realignunwarp.uweoptions.uwfwhm = 4;
    realignunwarp.uweoptions.rem = 1;
    realignunwarp.uweoptions.noi = 5;
    realignunwarp.uweoptions.expround = 'Average';
    realignunwarp.uwroptions.uwwhich = [2 1];
    realignunwarp.uwroptions.rinterp = 4;
    realignunwarp.uwroptions.wrap = WrapD;
    realignunwarp.uwroptions.mask = 1;
    realignunwarp.uwroptions.prefix = 'u';
    
    spm_run_realignunwarp(realignunwarp);

    rp_file = spm_file(funcfile, 'prefix','rp_','ext','.txt');

    keepfiles{numel(keepfiles)+1} = {spm_file(funcfile, 'prefix','rp_','ext','.txt')};

    reffunc = spm_file(funcfile, 'prefix','meanu');
    funcfile = spm_file(funcfile, 'prefix','u');

    delfiles{numel(delfiles)+1} = {reffunc};
    delfiles{numel(delfiles)+1} = {funcfile};
else
    %% Reslice the func series
    
    realignestwrite.data{1} = {funcfile};
    realignestwrite.eoptions.quality = 0.9;
    realignestwrite.eoptions.sep = 4;
    realignestwrite.eoptions.fwhm = 5;
    realignestwrite.eoptions.rtm = 1;
    realignestwrite.eoptions.interp = 2;
    realignestwrite.eoptions.wrap = [0 0 0];
    realignestwrite.eoptions.weight = '';
    realignestwrite.roptions.which = [2 1];
    realignestwrite.roptions.interp = 4;
    realignestwrite.roptions.wrap = [0 0 0];
    realignestwrite.roptions.mask = 1;
    realignestwrite.roptions.prefix = 'r';
    
    spm_run_realign(realignestwrite);

    rp_file = spm_file(funcfile, 'prefix','rp_','ext','.txt');

    keepfiles{numel(keepfiles)+1} = {spm_file(funcfile, 'prefix','rp_','ext','.txt')};

    reffunc = spm_file(funcfile, 'prefix','mean');
    funcfile = spm_file(funcfile, 'prefix','r');

    delfiles{numel(delfiles)+1} = {reffunc};
    delfiles{numel(delfiles)+1} = {funcfile};
end
function [rp_file,keepfiles] = my_spmbatch_acompcor(subfmridir,substring,task,wfuncdat,funcfile,rp_file,c1im,c2im,c3im,params,keepfiles)

fprintf('Start aCompCor \n')

if exist(rp_file,'file')
    confounds = load(rp_file);
else
    confounds = [];
end

GM = spm_vol(c1im);
WM = spm_vol(c2im);
CSF = spm_vol(c3im);

gmdat = spm_read_vols(GM);
wmdat = spm_read_vols(WM);
csfdat = spm_read_vols(CSF);

braindat = gmdat+wmdat;
braindat(braindat<0.2)=0;
braindat(braindat>0.0)=1;

csfdat(braindat>0.0)=0;
csfdat(csfdat<0.8)=0;
csfdat(csfdat>0.0)=1;

if params.do_bpfilter
    if params.nechoes==1
        funcjsonfile = fullfile(subfmridir,[substring '_task-' task '_bold.json']);
    else
        funcjsonfile = fullfile(subfmridir,[substring '_task-' task '_bold_e1.json']);
    end

    jsondat = fileread(funcjsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    bpfilter = [tr params.bpfilter(1:2)];
else
    bpfilter = [];
end

acc_confounds = fmri_acompcor(wfuncdat(:,:,:,:),{csfdat},params.Ncomponents,'confounds',confounds,'filter',bpfilter,'PolOrder',1);

if exist(rp_file)
    confounds = cat(2,confounds,acc_confounds);

    rp_file = spm_file(rp_file, 'prefix','acc_','ext','.txt');

    writematrix(confounds,rp_file,'Delimiter','tab');
else
    rp_file = spm_file(funcfile, 'prefix','acc_','ext','.txt');

    writematrix(acc_confounds,rp_file,'Delimiter','tab');
end

keepfiles{numel(keepfiles)+1} = {rp_file};

fprintf('Done aCompCor \n')
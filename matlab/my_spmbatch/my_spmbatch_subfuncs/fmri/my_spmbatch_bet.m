function [outfile,delfiles] = my_spmbatch_bet(infolder,infile,ppparams,params,delfiles,keepfiles)

%% Make GM, WM masks
if params.func.isaslbold
    if ~isfield(ppparams.perf(1),'c1m0scanfile') || ~isfield(ppparams.perf(1),'c2m0scanfile') || ~isfield(ppparams.perf(1),'c3m0scanfile')
        [ppparams,delfiles,keepfiles] = my_spmbatch_asl_segmentation(ppparams,params,delfiles,keepfiles);
    end

    c1im = fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile);
    c2im = fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile);
    c3im = fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile);
else
    if ~isfield(ppparams,'fc1im') || ~isfield(ppparams,'fc2im') || ~isfield(ppparams,'fc3im')
        fprintf('Do segmentation \n')
        
        preproc.channel.vols = {fullfile(ppparams.subfuncdir,reffunc)};
        preproc.channel.biasreg = 0.001;
        preproc.channel.biasfwhm = 60;
        preproc.channel.write = [0 0];
        preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')};
        preproc.tissue(1).ngaus = 1;
        preproc.tissue(1).native = [1 0];
        preproc.tissue(1).warped = [0 0];
        preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
        preproc.tissue(2).ngaus = 1;
        preproc.tissue(2).native = [1 0];
        preproc.tissue(2).warped = [0 0];
        preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
        preproc.tissue(3).ngaus = 2;
        preproc.tissue(3).native = [1 0];
        preproc.tissue(3).warped = [0 0];
        preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
        preproc.tissue(4).ngaus = 3;
        preproc.tissue(4).native = [0 0];
        preproc.tissue(4).warped = [0 0];
        preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
        preproc.tissue(5).ngaus = 4;
        preproc.tissue(5).native = [0 0];
        preproc.tissue(5).warped = [0 0];
        preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
        preproc.tissue(6).ngaus = 2;
        preproc.tissue(6).native = [0 0];
        preproc.tissue(6).warped = [0 0];
        preproc.warp.mrf = 1;
        preproc.warp.cleanup = 1;
        preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        preproc.warp.affreg = 'mni';
        preproc.warp.fwhm = 0;
        preproc.warp.samp = 3;
        preproc.warp.write = [0 1];
        preproc.warp.vox = NaN;
        preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];
        
        spm_preproc_run(preproc);
        
        ppparams.fc1im = fullfile(ppparams.subfuncdir,['c1' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        ppparams.fc2im = fullfile(ppparams.subfuncdir,['c2' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        ppparams.fc3im = fullfile(ppparams.subfuncdir,['c3' ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile]);
        
        sname = split(reffunc,'.nii');
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[sname{1} '._seg8.mat'])};
        delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['y_' sname{1} '.nii'])};  
    end
    
    c1im = ppparams.fc1im;
    c2im = ppparams.fc2im;
    c3im = ppparams.fc3im;
end

c1dat = spm_read_vols(spm_vol(c1im));
c2dat = spm_read_vols(spm_vol(c2im));
c3dat = spm_read_vols(spm_vol(c3im));

mask = c1dat+c2dat+c3dat;

Vin = spm_vol(fullfile(infolder,infile));
indat = spm_read_vols(Vin(1));

indat(find(mask<0.9)) = 0;

Vout = Vin(1);
Vout.fname = fullfile(infolder,['b' infile]);
Vout = spm_write_vol(Vout,indat);

outfile = ['b' infile];

delfiles{numel(delfiles)+1} = {fullfile(infolder,outfile)};

clear c1dat c2dat c3dat mask indat Vin Vout
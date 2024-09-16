function [delfiles,keepfiles] = my_spmbatch_preprocfasl(ppparams,params,delfiles,keepfiles)

%% Make GM, WM masks
if ~isfield(ppparams.perf(1),'c1m0scanfile') || ~isfield(ppparams.perf(1),'c2m0scanfile') || ~isfield(ppparams.perf(1),'c3m0scanfile')
    preproc.channel.vols = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile])};
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
    
    ppparams.perf(1).c1m0scanfile = ['c1' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];
    ppparams.perf(1).c2m0scanfile = ['c2' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];
    ppparams.perf(1).c3m0scanfile = ['c3' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile];

    fname = split([ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile],'.nii');
     
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[fname{1} '._seg8.mat'])};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,['y_' fname{1} '.nii'])};  
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile)};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile)};
end

%% Control-Label subtraction to make deltam series
if ~isfield(ppparams.perf(1),'deltamfile')
    fprintf('Start subtraction \n')

    [ppparams,delfiles,keepfiles] = my_spmbacth_faslsubtraction(ppparams,params,delfiles,keepfiles);
end

%% Make cbf series
if ~isfield(ppparams.perf(1),'cbffile')
    fprintf('Start CBF mapping \n')
    
    [ppparams,delfiles,keepfiles] = my_spmbatch_fasl_cbfmapping(ppparams,params,delfiles,keepfiles);
end

%% Normalise CBF series
if ~isfield(ppparams.perf(1),'wcbffile')
    fprintf('Do normalization\n')

    [ppparams,delfiles,keepfiles] = my_spmbatch_aslbold_normalization(ppparams,params,delfiles,keepfiles);
end

%% Smooth CBF series        
if ~isfield(ppparams.perf(1),'scbffile')
    fprintf('Do smoothing \n')

    [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothasl(ppparams,params,delfiles,keepfiles);
end
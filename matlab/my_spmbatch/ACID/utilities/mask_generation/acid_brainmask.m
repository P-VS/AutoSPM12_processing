function Ni = acid_brainmask(P, P_others, ~, dummy_options, dummy_segmentation)

% =========================================================================
% This function maskes any image PBMSK. If the dummy_options are used can be used for each kind of image (mean DW or b0)
% Created by S.Mohammadi 08.07.2015
% =========================================================================

V = spm_vol(P);

% define output directory
keyword = 'BRAIN-MASK';
[path,fname,~] = spm_fileparts(V.fname);

p_out = acid_bids(path,fname,keyword,1);


if dummy_segmentation == 1

    I_others = acid_read_vols(V,V,1);
    
    P = acid_write_vol(I_others, V, p_out, 'BRAIN-TPM');
    P = P.fname;

    [spm_path,~,~] = fileparts(which('spm.m'));
    tpm_path = [spm_path filesep 'tpm' filesep 'TPM.nii'];

    % create TPMs
    segment_job.channel.vols = {P};
    segment_job.channel.biasreg = 0.001;
    segment_job.channel.biasfwhm = 60;
    segment_job.channel.write = [0 0];
    segment_job.tissue(1).tpm = {[tpm_path ',1']};
    segment_job.tissue(1).ngaus = 1;
    segment_job.tissue(1).native = [1 0];
    segment_job.tissue(1).warped = [0 0];
    segment_job.tissue(2).tpm = {[tpm_path ',2']};
    segment_job.tissue(2).ngaus = 1;
    segment_job.tissue(2).native = [1 0];
    segment_job.tissue(2).warped = [0 0];
    segment_job.tissue(3).tpm = {[tpm_path ',3']};
    segment_job.tissue(3).ngaus = 2;
    segment_job.tissue(3).native = [1 0];
    segment_job.tissue(3).warped = [0 0];
    segment_job.tissue(4).tpm = {[tpm_path ',4']};
    segment_job.tissue(4).ngaus = 3;
    segment_job.tissue(4).native = [1 0];
    segment_job.tissue(4).warped = [0 0];
    segment_job.tissue(5).tpm = {[tpm_path ',5']};
    segment_job.tissue(5).ngaus = 4;
    segment_job.tissue(5).native = [1 0];
    segment_job.tissue(5).warped = [0 0];
    segment_job.tissue(6).tpm = {[tpm_path ',6']};
    segment_job.tissue(6).ngaus = 2;
    segment_job.tissue(6).native = [0 0];
    segment_job.tissue(6).warped = [0 0];
    segment_job.warp.mrf = 1;
    segment_job.warp.cleanup = 1;
    segment_job.warp.reg = [0 0.001 0.5 0.05 0.2];
    segment_job.warp.affreg = 'mni';
    segment_job.warp.fwhm = 0;
    segment_job.warp.samp = 3;
    segment_job.warp.write = [0 0];
    segment_job.warp.vox = NaN;
    segment_job.warp.bb = [NaN NaN NaN NaN NaN NaN];

    segment_struct = spm_preproc_run(segment_job);

    % rename TPMs
    cd(p_out);
    [~,c1,~] = spm_fileparts(char(segment_struct(1).tiss(1).c));
    [~,c2,~] = spm_fileparts(char(segment_struct(1).tiss(2).c));
    [~,c3,~] = spm_fileparts(char(segment_struct(1).tiss(3).c));
    [~,c4,~] = spm_fileparts(char(segment_struct(1).tiss(4).c));
    [~,c5,~] = spm_fileparts(char(segment_struct(1).tiss(5).c));

    c1_new = [c1(3:end) '1.nii'];
    c2_new = [c2(3:end) '2.nii'];
    c3_new = [c3(3:end) '3.nii'];
    c4_new = [c4(3:end) '4.nii'];
    c5_new = [c5(3:end) '5.nii'];

    movefile([c1 '.nii'], c1_new);
    movefile([c2 '.nii'], c2_new);
    movefile([c3 '.nii'], c3_new);
    movefile([c4 '.nii'], c4_new);
    movefile([c5 '.nii'], c5_new);

    P_tpm(1,:) = [p_out filesep c1_new];
    P_tpm(2,:) = [p_out filesep c2_new];
    P_tpm(3,:) = [p_out filesep c3_new];

    delete([c1_new(1:end-5) '_seg8.mat']);
    clear P_mask
end

if dummy_segmentation == 1
    V_tpm(1,:) = spm_vol([P_tpm(1,:) ',1']);
    V_tpm(2,:) = spm_vol([P_tpm(2,:) ',1']);
    V_tpm(3,:) = spm_vol([P_tpm(3,:) ',1']);
else
    V_tpm = spm_vol(P_tpm);
end

% This function smoothes first the segmented image to ensure that only
% voxels within the tissue segments are included. Furthermore, it
% automatically increases the threshold to cover the brain.
I_tpm = sum(spm_read_vols(V_tpm),4);
smk   = dummy_options.smk;
perc  = dummy_options.perc;
deltainx = 0.01;

% smooth segments
spm_smooth(I_tpm,I_tpm,smk);

% determine threshold for mask
[y,x] = hist(I_tpm(:),100);
cy    = cumsum(y);
sz    = size(I_tpm(:),1);
THR   = x(max(find(cy<=sz*perc)));

while isempty(THR)
    if(perc>=0.99)
        deltainx = 1e-3;
    end
    perc = perc+deltainx;
    THR  = x(max(find(cy<=sz*perc)));
    disp(['Threshold for brain mask: ' num2str(perc)])
end

% generate filename
V_mask  = V_tpm(1);
V_mask.fname = V.fname;
fname = acid_bids_filename(V_mask, keyword, 'same', '.nii');
fname = [p_out filesep fname];

% write out binary brain mask
Ni = nifti;
Ni.mat  = V_mask(1).mat;
Ni.mat0 = V_mask(1).mat;
vol  = zeros(V_mask.dim);
vol(find(I_tpm>THR)) = 1;
dt  = [spm_type('float32'),spm_platform('bigend')];
Ni.dat = file_array(fname,size(vol),dt, 0,1,0);
create(Ni);

spm_progress_bar('Init',size(vol,4),Ni.descrip,'volumes completed');

for p = 1:size(vol,4)
    Ni.dat(:,:,:,p) = vol(:,:,:,p);
    spm_progress_bar('Set',p);
end

spm_progress_bar('Clear');

if ~isempty(P_others)
    acid_brainmask_apply(fname, P_others, p_out)
end

delete([c1_new(1:end-5) '.nii']);

end
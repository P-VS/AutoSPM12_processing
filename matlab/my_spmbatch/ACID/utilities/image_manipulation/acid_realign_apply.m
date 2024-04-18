function V_out = acid_realign_apply(PT, PS, params)

% =========================================================================
%
% Inputs:
%   PT     - string containing path and name of target image
%   PS     - string containing path and name of other images to apply slicewise affine transformation
%   params - array containing transformation parameters for each slide
%   p_out  - output directory
%
% Outputs:
%   V_out  - header of the resliced image
%
% Created by S.Mohammadi 12/02/2016
% =========================================================================

% interpolation kernel
interp_def = -4;

for i = 1:size(params,1)
    load(deblank(params(i,:)));
    paramsslice{i} = params;
end

% load in target image
VT = spm_vol(PT);
dm = VT.dim;

% load in source image(s)
if size(PS,1)
    VS = spm_vol(PS);
else
    VS = acid_load_4Dimage(PS);
end

% define output directory
keyword = 'REALIGN';
[path,fname,~] = spm_fileparts(VS(1).fname);

p_out = acid_bids(path,fname,keyword,1);


% save json file
acid_save_json(VS, p_out, keyword);

% save bvals and bvecs files
fname = acid_bids_filename(VS,keyword,'_dwi','.bval');
fname = [p_out filesep fname];
dlmwrite(fname, bvals_new);

fname = acid_bids_filename(VS,keyword,'_dwi','.bvec');
fname = [p_out filesep fname];
dlmwrite(fname, bvecs_new);

    % save json file
    acid_save_json(VO(1), p_out, keyword);
    
    % save bvals and bvecs files
    acid_save_bvals_bvecs(VO(1), p_out, keyword);

% transform parameter file
paramout = zeros(size(params));
for z = 1:dm(3)
    mat = eye(4);
    for i = 1:numel(paramsslice)
        paramtmp = paramsslice{i};
        mat = spm_matrix(paramtmp(z,:))*mat;
    end
    paramout(z,:) = spm_imatrix(mat);
end

% apply transformation to source images
sz = size(VS,1);
for inx = 1:sz
    I = zeros(dm);
    V_out = VT;
    V_out.fname = VS(inx).fname;
    for p=1:dm(3)
        xpar    = paramout(p,:);
        iM      = inv(spm_matrix_mod(xpar));
        M       = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(VT.mat)*iM*VS(inx).mat);
        tmp     = spm_slice_vol(VS(inx),M,dm(1:2),interp_def);
        I(:,:,p) = tmp;
    end
    I_out(:,:,:,inx) = I;
end

    % write realigned images
    V_out = acid_write_vol(I_out, Vout, p_out, keyword, 'same');
end
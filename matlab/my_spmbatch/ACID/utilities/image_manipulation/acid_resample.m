function VO = acid_resample(P, res, BB, dummy_mask, interp)

% ========================================================================
% acid_resample -- resample images to have specified voxel dims and BBox
% acid_resample(imnames, res, bb, dummy_mask)
%
% Output images will be prefixed with 'r', and will have voxel dimensions
% equal to res. Use NaNs to determine ress from transformation matrix
% of input image(s).
% If bb == nan(2,3), bounding box will include entire original image
% Origin will move appropriately. Use world_bb to compute bounding box from
% a different image.
%
% Pass dummy_mask=true to re-round binary mask values (avoid
% growing/shrinking masks due to linear interp)
%
% See also res, world_bb
%
% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
% Adapted by Ged Ridgway -- email bugs to drc.spm@gmail.com

% This version doesn't check spm_flip_analyze_images -- the handedness of
% the output image and matrix should match those of the input.
%
% modified by S.Mohammadi 06/05/2014
% modified by B.Fricke 07/03/2022
% To ensure that the same voxel-sign is used after interpolation
% To distinguish this step from coreg, if is written out with an "i"
% ========================================================================

% Check spm version:
if exist('spm_select','file') % should be true for spm5
    spm5 = 1;
elseif exist('spm_get','file') % should be true for spm2
    spm5 = 0;
else
    error('Can''t find spm_get or spm_select; please add SPM to path')
end

spm_get_defaults;

% prompt for missing arguments
if ~exist('P','var') || isempty(char(P))
    if spm5
        P = spm_select(inf, 'image', 'Choose images to resize');
    else
        P = spm_get(inf, 'img', 'Choose images to resize');
    end
end

% check if inter fig already open, don't close later if so...
Fint = spm_figure('FindWin', 'Interactive'); Fnew = [];
if ~exist('res', 'var') || isempty(res)
    Fnew = spm_figure('GetWin', 'Interactive');
    res = spm_input('Vox Dims (NaN for "as input")? ',...
        '+1', 'e', '[nan nan nan]', 3);
end
if ~exist('BB', 'var') || isempty(BB)
    Fnew = spm_figure('GetWin', 'Interactive');
    BB = spm_input('Bound Box (NaN => original)? ',...
        '+1', 'e', '[nan nan nan; nan nan nan]', [2 3]);
end
if ~exist('dummy_mask', 'var')
    dummy_mask = false;
end
if isempty(dummy_mask)
    dummy_mask = false;
end
if ~exist('interp', 'var')
    % Default to 7th order sinc interpolation:
    interp = -7;
end
if isempty(interp)
    interp = -7;
end


% return error for a single(3D) input
if size(P,1) == 1
    struct4D = nifti(P);
    dim4D = struct4D.dat.dim;
    try n = dim4D(4); 
    catch 
    end
    if length(dim4D) < 4
            V = spm_vol(P);
    elseif n == 1
            V = spm_vol(P);
    else
            % load 4D image
            V = acid_load_4Dimage(P);
    end
end


% define output directory
keyword = 'RESAMP';
[path,fname,~] = spm_fileparts(V(1).fname);

p_out = acid_bids(path,fname,keyword,1);


j = 1;

for V = V'
    
    % copy to allow defaulting of NaNs differently for each volume
    res = sign(ones(1,3)*V.mat(1:3,1:3)).*abs(res); % SM
    bb = BB;
    
    % default res to current volume's res, (from mat parameters)
    if any(isnan(res))
        vprm = spm_imatrix(V.mat);
        vres = vprm(7:9);
        res(isnan(res)) = vres(isnan(res));
    end
    res = res(:)';

    mn = bb(1,:);
    mx = bb(2,:);
    % default BB to current volume's
    if any(isnan(bb(:)))
        vbb = world_bb(V);
        vmn = vbb(1,:);
        vmx = vbb(2,:);
        % SM 
        sgn = sign(res);
        vpos = find(sgn == -1);
        for inx = vpos
            vmn(inx) = vbb(2,inx);
            vmx(inx) = vbb(1,inx);
        end
         % SM
        mn(isnan(mn)) = vmn(isnan(mn));
        mx(isnan(mx)) = vmx(isnan(mx));
    end

    % voxel [1 1 1] of output should map to BB mn
    % (the combination of matrices below first maps [1 1 1] to [0 0 0])
    mat = spm_matrix([mn 0 0 0 res])*spm_matrix([-1 -1 -1]);
    % voxel-coords of BB mx gives number of voxels required
    % (round up if more than a tenth of a voxel over)
    
    % SM
    % imgdim = ceil(mat \ [mx 1]' - 0.1)';
    imgdim = abs(ceil(mat \ [mx 1]' - 0.1)'); % SM
    imgdim = imgdim - 1;
    % SM
    
    % output image
    VO     = V;
    VO.mat = mat;
    spm_progress_bar('Init',imgdim(3),'reslicing...','planes completed');
    for i = 1:imgdim(3)
        M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
        img = spm_slice_vol(V, M, imgdim(1:2), interp);
        if dummy_mask
            img = round(img);
        end
        V_out(:,:,i) = img;
        spm_progress_bar('Set', i)
    end
    
    I_out(:,:,:,j) = V_out;
    j = j+1;
    spm_progress_bar('Clear');
end

% write resampled image
VO = acid_write_vol(I_out, VO, p_out, keyword, 'same', '', '', 1);

% save json file
acid_save_json(V(1), p_out, keyword);

% save bvals and bvecs
acid_save_bvals_bvecs(V, p_out, keyword);

% closing interactive figure opened by this script
if (isempty(Fint) && ~isempty(Fnew)) 
    close(Fnew);
end

disp('Done.')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bb = world_bb(V)
    %  world-bb -- get bounding box in world (mm) coordinates

    d = V.dim(1:3);
    % corners in voxel-space
    c = [ 1    1    1    1
          1    1    d(3) 1
          1    d(2) 1    1
          1    d(2) d(3) 1
          d(1) 1    1    1
          d(1) 1    d(3) 1
          d(1) d(2) 1    1
          d(1) d(2) d(3) 1 ]';
    % corners in world-space
    tc = V.mat(1:3,1:4)*c;

    % bounding box (world) min and max
    mn = min(tc,[],2)';
    mx = max(tc,[],2)';
    bb = [mn; mx];
end
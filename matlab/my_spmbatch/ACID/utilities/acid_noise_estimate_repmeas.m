function [sigma,SNRout,masknames] = acid_noise_estimate_repmeas(P, P_mask, bvals, dummy_shell, rounding_tolerance)

%==========================================================================
% Estimates SNR from repeated measures using the formula in Dietrich et al., JMRI, 2007;  DOI: 10.1002/jmri.20969
%
% Inputs:
% P      - file names of DTI images (i.e. low and high b-value images 
% P_mask - region of interest, in which the tensor estimates are calculated 
% bvals  - b-values of the source images
%
% Created by S.Mohammadi 31/05/2021
%==========================================================================
tic

% prepare if first volume of a 4D file was selected
if size(P,1) == 1
    struct4D = nifti(P);
    dim4D = struct4D.dat.dim;
    try n = dim4D(4);
       if n == 1
           error('A single 3D source image was selected. Choose a single 4D volume or select all 3D volumes manually!');
       end
       P = strcat(repmat(P(1:end-2), n, 1), ',', num2str((1:n)'));
       catch
   end
end

V      = spm_vol(P);
V_mask = spm_vol(P_mask);
SNRout = zeros(1,size(V_mask,1));
sigma  = zeros(1,size(V_mask,1));
masknames = cell(1,size(V_mask,1));

keyword = 'SIGMA';
[path,fname,~] = spm_fileparts(V(1).fname);

p_out = acid_bids(path,fname,keyword,1);
diary([p_out filesep 'logfile_' keyword '.txt'])

% defaults
res = -3;

bvals = round(bvals/rounding_tolerance)*rounding_tolerance;

b0     = min(bvals);
MSK_b0 = find(bvals == b0);
MSK_b  = find(bvals > b0);

if size(MSK_b0,2) > (size(bvals,2)/2)
    warning('More than 50% of your b-values were set to zero! Plesase check whether you used um/ms as units and did not adjust the rounding_tolerance in the defaults!');
end

[b_unique,~,~] = unique(bvals(MSK_b));

% n = histc(bvals(MSK_b),b_unique);

if dummy_shell == 1
    selected_shell = max(b_unique);
else
    selected_shell = 0;
end

disp(['Measured shells:                    ' num2str(b_unique)]);
disp(['Selected shell for noise estimation: b' num2str(selected_shell)]);


% make bval mask
mskbval = find(bvals==selected_shell);

if size(mskbval,2) < 2
    error('At least two volumes of the same shell (or b0) are needed for the repeated measurement methode.');
end

% function SNR of an ROI for multiple measurement
fSNR   = @(yIn) mean(yIn(:))/mean(std(yIn,[],1));
fnoise = @(yIn) mean(std(yIn,[],1));



for inx_msk = 1:size(V_mask,1)
    Amsk    = acid_read_vols(V_mask(inx_msk),V(1),res);
    msk     = find(Amsk>0);
    yIn     = zeros(length(mskbval),numel(msk));
    
    for inx = 1:length(mskbval)
        A = acid_read_vols(V(mskbval(inx)),V(1),res);
        yIn(inx,:) = A(Amsk>0);
    end
    
    SNRout(inx_msk) = fSNR(yIn);
    sigma(inx_msk)  = fnoise(yIn);
    [~,f,~] = spm_fileparts(V_mask(inx_msk).fname);
    masknames{inx_msk} = f;
end

% save sigma estimates

fname = acid_bids_filename(V(1),keyword);
fname_sigma = [p_out filesep fname(1:end-4) '.txt'];
save(fname_sigma,'sigma','-ascii');

disp(['Sigma: ' num2str(sigma)])
T = toc/60;

T = duration(minutes(T),'format','hh:mm:ss.SSS');
disp(['The total time for ' keyword ' was: ' char(T) '.']);
diary off

end
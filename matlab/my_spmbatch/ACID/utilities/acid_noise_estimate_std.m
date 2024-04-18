function [sigma, sigma_avg] = acid_noise_estimate_std(P, P_mask, ncoil)

%==========================================================================
% Estimates SNR from noise mask using the formula in Hutton et al., Plos ONE, 2012 (doi:10.1371/journal.pone.0052075).
%
% Inputs:
% P      - file names of DTI images (i.e. low and high b-value images 
% P_mask - region of interest, in which the tensor estimates are calculated 
% ncoil  - effective number of coils
%
% Created by S.Mohammadi 11/09/2013
%==========================================================================
tic
if ~exist('P','var')
    P = char(cfg_getfile(Inf,'IMAGE','Get images','','','.*'));
end

if ~exist('P_mask','var')
    P_mask = char(cfg_getfile(1,'IMAGE','Get ROI mask','','','.*'));
end

if ~exist('ncoil','var')
    ncoil = spm_input('Select number of coils',1);
end

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

% load in dwi data
V = spm_vol(P);

keyword = 'SIGMA';
[path,fname,~] = spm_fileparts(V(1).fname);

p_out = acid_bids(path,fname,keyword,1);
diary([p_out filesep 'logfile_' keyword '.txt'])
% load in mask and check dimension
if ~isempty(P_mask)
    V_mask = spm_vol(P_mask);
    if(~isequal(V_mask.dim,V(1).dim))
        error('Dimension of Mask and DTIs must agree!')
    end
    I_mask = spm_read_vols(V_mask); 
end

% estimate sigma
sigma = zeros(length(V),1);
for i = 1:size(P,1)
    I = spm_read_vols(V(i));
    I_maskout    = (I_mask>0 & ~isnan(I));
    sigma(i) = sqrt(mean(I(I_maskout).^2/2/ncoil));
end

sigma_avg = mean(sigma);

% save sigma estimates


fname = acid_bids_filename(V(1),keyword);
fname_sigma = [p_out filesep fname(1:end-4) '.txt'];
save(fname_sigma,'sigma','-ascii');

fname_sigma_avg = [p_out filesep fname(1:end-4) '-avg' '.txt'];
save(fname_sigma_avg,'sigma_avg','-ascii');

disp(['Sigma: ' num2str(sigma_avg)])
T = toc/60;

T = duration(minutes(T),'format','hh:mm:ss.SSS');
disp(['The total time for ' keyword ' was: ' char(T) '.']);
diary off
end
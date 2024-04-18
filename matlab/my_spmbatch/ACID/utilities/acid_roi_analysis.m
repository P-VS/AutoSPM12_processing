function ROI_avg = acid_roi_analysis(P_maps, P_maskGLO, P_maskSUB, P_maskREL, dummy_slicewise)

% ===========================================================================
% The following script extracts average values within masks
%
% Inputs:
%   P_maps  - input maps
%
%   P_masks   - subject-specific masks
%                 -> requires a character array with filenames of the masks;
%                    if not specified, file selector pops up
%   P_maskGLO - global mask
%                 -> requires the filename of the global mask
%                    if not specified, file selector pops up
%
%   P_maskSUB - subject specific masks
%                 -> requires a character array with filenames of the masks;
%                    if not specified, file selector pops up
%
%   P_maskREL - subject specific reliability masks
%                 -> requires a character array with filenames of the reliability masks;
%                    if not specified, file selector pops up
%
%   P_dummy_slicewise - specifies whether the script extract values from the
%                   whole volume or in each slice separately
%                   0: volume-wise extraction
%                   1: slice-wise extraction
%
%   p_out - output directory
%
% Output:
%   - avgROI: average value within the intersection of masks specified above
%
% Created by G.David
% ===========================================================================

% check for existence of input parameters
if ~exist('P_maps','var')
    P_maps = char(cfg_getfile([1 Inf],'any','Select DTI maps','','','.*'));
end
if ~exist('P_maskGLO','var')
    P_maskGLO = char(cfg_getfile([1 Inf],'any','Select ROIs','','','.*'));
end
if isempty(P_maskGLO)
    dummy_maskGLO = 0;
else
    dummy_maskGLO = 1;
end
if ~exist('P_maskSUB','var')
    P_maskSUB = char(cfg_getfile([0 Inf],'any','Select subject-specific masks','','','.*'));
end
if isempty(P_maskSUB)
    dummy_maskSUB = 0;
else
    dummy_maskSUB = 1;
    if size(P_maskSUB,1)~=size(P_maps,1)
        error('The number of specified subject-specific ROIs has to match the number of DTI maps.');
    end
end
if ~exist('P_maskREL','var')
    P_maskREL = char(cfg_getfile([0 Inf],'any','Select reliability masks','','','.*'));
end
if isempty(P_maskREL)
    dummy_maskREL = 0;
else
    dummy_maskREL = 1;
    if size(P_maskREL,1)~=size(P_maps,1)
        error('The number of specified reliability masks has to match the number of DTI maps.');
    end
end
if ~exist('dummy_slicewise','var')
    dummy_slicewise = spm_input('Slice-wise average? 0 (no), 1 (yes)',1);
end

% initialization
if ~dummy_slicewise
    if ~dummy_maskGLO
        ROI_avg = nan(size(P_maps,1),1);
        ROI_std = nan(size(P_maps,1),1);
    else
        ROI_avg = nan(size(P_maps,1),size(P_maskGLO,1));
        ROI_std = nan(size(P_maps,1),size(P_maskGLO,1));
    end
else
    sl = zeros(size(P_maps,1),1);
    for  i = 1:size(P_maps,1)
        V_map = spm_vol(P_maps(i,:));
        sl(i) = V_map.dim(3);
    end
    
    if ~dummy_maskGLO
        ROI_avg = nan(size(P_maps,1),max(sl));
        ROI_std = nan(size(P_maps,1),max(sl));
    else
        ROI_avg = nan(size(P_maps,1),size(P_maskGLO,1),max(sl));
        ROI_std = nan(size(P_maps,1),size(P_maskGLO,1),max(sl));
    end
end

for k = 1:size(P_maps,1)
    
    disp(['Processing: ' P_maps(k,:)]);
    
    % load in maps
    V_map = spm_vol(P_maps(k,:));
    maps = acid_read_vols(V_map,V_map,1);
    dm = V_map.dim;
           
    % masking with subject-specific ROI
    if dummy_maskSUB
        V_maskSUB = spm_vol(P_maskSUB(k,:));
        MaskSUB = acid_read_vols(V_maskSUB,V_map,0);
        MaskSUB(MaskSUB>10^-3) = 1;
        maps = maps.*MaskSUB;
    end
    
    % masking with subject-specific reliability mask
    if dummy_maskREL
        V_maskREL = spm_vol(P_maskREL(k,:));
        MaskREL = acid_read_vols(V_maskREL,V_map,0);
        MaskREL(MaskREL>10^-3) = 1;
        maps = maps.*MaskREL;
    end
    
    % global ROIs
    if dummy_maskGLO
        for m = 1:size(P_maskGLO,1)
            V_maskGLO = spm_vol(P_maskGLO(m,:));
            ROI = acid_read_vols(V_maskGLO,V_map,0);
            ROI(ROI>10^-3)=1;
            maps2 = maps.*ROI;
            if ~dummy_slicewise
                temp = maps2; temp = temp(:); temp(temp==0) = []; temp(isnan(temp)) = [];
                ROI_avg(k,m) = mean(temp);
                ROI_std(k,m) = std(temp);
            else
                for i = 1:dm(3)
                    temp = squeeze(maps2(:,:,i));  temp = temp(:); temp(temp==0) = []; temp(isnan(temp)) = [];
                    ROI_avg(k,m,i) = mean(temp);
                    ROI_std(k,m,i) = std(temp);
                end
            end
        end
    else
        if ~dummy_slicewise
            temp = maps; temp = temp(:); temp(temp==0) = []; temp(isnan(temp)) = [];
            ROI_avg(k,1) = mean(temp);
            ROI_std(k,1) = std(temp);
        else
            for i = 1:dm(3)
                temp = squeeze(maps(:,:,i));  temp = temp(:); temp(temp==0) = []; temp(isnan(temp)) = [];
                ROI_avg(k,i) = mean(temp);
                ROI_std(k,i) = std(temp);
            end
        end
    end    
end

% define output directory
[p_out,~,~] = fileparts(P_maps(1,:));

fout_mean = [p_out filesep 'ROI_mean.mat'];
fout_std = [p_out filesep 'ROI_std.mat'];

save(fout_mean, 'ROI_avg');
save(fout_std,  'ROI_std');

disp(['ROI average: ' num2str(ROI_avg)]);
disp(['ROI standard deviation: ' num2str(ROI_std)]);

end
function thr = acid_relmask_thr(errmaps, maps, masks, maskGLO, emin, emax, nstep, dummy_display)

% ================================================================================
% The function determines the optimal threshold for reliability masking.
% The procedure is described in detail in David et al., 2017.
% In short, the method is based on minimizing the standard error of the
% mean (sem) of FA across a homogeneous pool of voxels.
%
% Inputs:
%   - errmaps: model-fit error maps corresponding to the FA maps specified above;
%                the script expects a char array with the filenames of
%                model-fit error maps in the rows;
%                (hint: use spm_select to generate this input)
%                Note: model-fit error here refers to the rms(error) within a given voxel
%
%   - maps: a set of FA (or other DTI scalar) maps;
%           the script expects a char array with the filenames of the FA
%           maps in the rows;
%           (hint: use spm_select to generate this input)
%
%   - masks: subject-specific masks
%               -> requires a character array with filenames of the masks;
%                  (hint: use spm_select to generate this input)
%
%   - maskGLO: global mask applied for all subjects
%               -> requires the filename of the global mask
%
%   - emin: emin*median(error) is the the lowest tested model-fit error for
%           reliability masking (default: 1)
%
%   - emax: emax*median(error) is the the highest tested model-fit error for
%           reliability masking (default: 3)
%
%   - nstep: how many equidistal values are tested between emin*median(error)
%            and emax*median(error) for reliability masking
%
%   - dummy_display: flag whether
%                   0: no plot
%                   1: plot std, number of voxels, sem
%
%   - p_out: output folder
%
% Outputs:
%   - thr: optimal threshold for reliability masking
%
% Created by G.David, May 2017
% ================================================================================

% check for compulsory input
if ~exist('maps','var')
    error('No DTI maps specified.');
elseif ~exist('errmaps','var')
    error('No model-fit error maps specified.');
end

% define dummies
if ~isempty(maskGLO)
    dummy_maskGLO = 1;
else
    dummy_maskGLO = 0;
end
if ~isempty(masks)
    dummy_masks = 1;
else
    dummy_masks = 0;
end
if ~exist('dummy_display','var')
    dummy_display = 0;
end

% optional input
if ~exist('emin','var')
    emin = 1;
end
if ~exist('emax','var')
    emax = 3;
end
if ~exist('nstep','var')
    nstep = 101;
end

if size(errmaps,1)~=size(maps,1)
    error('The number of model-fit error and FA maps has to be the same.');
elseif size(errmaps,1)~=size(masks,1)
    error('The number of model-fit error maps and subject-specific masks has to be the same.');
end

% load in reference map
VG = spm_vol(maps(1,:));

% load in global mask
if dummy_maskGLO
    I_maskGLO = spm_vol(maskGLO(1,:));
    maskGLO = acid_read_vols(I_maskGLO,I_maskGLO,1);   
end

% initialize stuff
map_group = cell(size(maps,1),1);
err_group = cell(size(maps,1),1);

for k = 1:size(maps,1)

    [folder,~,~] = fileparts(errmaps(k,:));
    fprintf('Processing subject: %s\n',folder);
    
    % load in DTI and corresponding model-fit error maps
    I_map = spm_vol(maps(k,:));
    I_err = spm_vol(errmaps(k,:));
    map = acid_read_vols(I_map,I_map,1);
    err = acid_read_vols(I_err,I_err,1);
    if ~isequal(size(err),size(map))
        error('The dimensions of the model-fit error and corresponding FA map have to be the same.');
    end
    
    % masking with subject-specific mask
    if dummy_masks
        I_mask = spm_vol(masks(k,:));
        mask = acid_read_vols(I_mask,I_mask,1);
        if ~isequal(size(err),size(mask))
            error('The dimensions of the model-fit error map and corresponding subject-specific mask have to be the same.');
        end
        map  = map.*mask;
        err  = err.*mask;
    end
    
    % masking with global mask
    if dummy_maskGLO
        map = map.*maskGLO;
        err = err.*maskGLO;
    end
    
    % vectorize
    map = map(:);
    err = err(:);
       
    % exclude voxels with 0 or NaN
    temp = map;
    map(temp==0)=[];
    err(temp==0)=[];
    
    msk_idx = (isnan(map) | isnan(err));
    map(msk_idx)=[];
    err(msk_idx)=[];
    
    % store remaining voxels in cell arrays
    map_group{k,1} = map;
    err_group{k,1} = err;
        
end

% pool across subjects
map_dist = map_group{1,1};
err_dist = err_group{1,1};
for k = 2:size(maps,1)
    map_dist = cat(1,map_dist,map_group{k,1});
    err_dist = cat(1,err_dist,err_group{k,1}); 
end

% determine median model-fit error
thr_median = median(err_dist);
thr_vect = linspace(emin,emax,nstep)*thr_median;

Ratio_n = zeros(length(thr_vect),1);
Ratio_std = zeros(length(thr_vect),1);
Ratio_sem = zeros(length(thr_vect),1);

% looping through thresholds
for i = 1:length(thr_vect)
    % mask image with REL mask.
    Mmap_dist = map_dist.*double(err_dist < thr_vect(i)); 
    Mmap_dist(Mmap_dist==0)=[];
    
    Ratio_n(i)    = (numel(Mmap_dist)-numel(map_dist))/numel(map_dist)*100;
    Ratio_std(i)  = (std(Mmap_dist)-std(map_dist))/std(map_dist)*100;
    Ratio_sem(i)  = (std(Mmap_dist)./sqrt(numel(Mmap_dist)) - std(map_dist)/sqrt(numel(map_dist))) / (std(map_dist)/sqrt(numel(map_dist)))*100;    
    Ratio_mean(i) = (mean(Mmap_dist) - mean(map_dist)) / mean(map_dist) *100;
    Ratio_t(i)    = (mean(Mmap_dist)./std(Mmap_dist).*sqrt(numel(Mmap_dist)) - mean(map_dist)/std(map_dist)*sqrt(numel(map_dist)) ) / ( mean(map_dist)/std(map_dist)*sqrt(numel(map_dist)) ) *100;
    
end

% optimal threshold
tmp = find(Ratio_sem==min(Ratio_sem));
thr = (1+(tmp(1)-1)*(emax-emin)/(nstep-1))*thr_median;

% define output directory
    [path,~,~] = spm_fileparts(VG.fname);
    p_out = path;


% save threshold reliability masking
save([p_out filesep 'threshold_reliability_masking.mat'],'thr');

% display results
if dummy_display
    figure
    plot(linspace(1,3,nstep),Ratio_std,'k',linspace(1,3,nstep),Ratio_n,'b',linspace(1,3,nstep),Ratio_sem,'r');
   %plot(linspace(1,3,nstep),Ratio_std,'k',linspace(1,3,nstep),Ratio_n,'b',linspace(1,3,nstep),Ratio_sem,'r',linspace(1,3,nstep),Ratio_mean,'g',linspace(1,3,nstep),Ratio_t,'y');
    temp = ylim;
    hold on
    plot(1+(tmp(1)-1)*(emax-emin)/(nstep-1),temp(1):temp(2),'k')
    title('determining optimal threshold')
    xlabel('rms model-fit error [multiples of median]')
    ylabel('relative change [%]')
    yline(0,'k')
    legend('std','sample size','sem')
    saveas(gcf, [p_out filesep 'threshold_reliability_masking.fig'])
end
end
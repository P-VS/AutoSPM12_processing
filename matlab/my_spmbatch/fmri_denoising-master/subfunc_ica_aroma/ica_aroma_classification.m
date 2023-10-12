function noiseICdata = ica_aroma_classification(ppparams, funcmask, ica_dir, t_r)

%This function classifies the ica components extracted by fmri_do_ica into
% noise/no-noise (both physiological and motion noise) based on:
%- The maximum robust correlation with the noise regressors (confounds):
%         - aCompCor confounds (physiological noise components) identified by aCompCor (the temporal signals of 5PCAcomponents determined
%           in the CSF using aCompCor, as physiological noise regressors.)
%         - 24 the motion regressors computed from the realignment
%           parameters
%- The amount of high frequency (> 1 Hz)
%- Their location at the brain edge

aroma_dir = fullfile(ica_dir,'ICA-AROMA');
if ~exist(aroma_dir,"dir"), mkdir(aroma_dir); end

GM = ppparams.wc1im;
WM = ppparams.wc2im;
CSF = ppparams.wc3im;

gmdat = spm_read_vols(spm_vol(GM));
wmdat = spm_read_vols(spm_vol(WM));
csfdat = spm_read_vols(spm_vol(CSF));

braindat = spm_read_vols(spm_vol(funcmask));
braindat((gmdat + wmdat) < 0.2) = 0;
braindat(braindat > 0.0) = 1;
    
Vbr = spm_vol(funcmask);
Vbr.fname = fullfile(aroma_dir, 'braindata.nii'); 
Vbr.descrip = 'braindata';
Vbr = rmfield(Vbr, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
spm_write_vol(Vbr, braindat);

nobrain = spm_read_vols(spm_vol(funcmask)); 
nobrain(braindat > 0.0) = 0; 
nobrain(nobrain > 0.0) = 1;

Vcsf = spm_vol(funcmask);
Vcsf.fname = fullfile(aroma_dir, 'nobraindata.nii'); 
Vcsf.descrip = 'nobraindata';
Vcsf = rmfield(Vcsf, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
spm_write_vol(Vcsf, nobrain);

edgedat = edge3(braindat,"approxcanny",0.2);
edgethn = ones(2,2,2);
edgedat = convn(edgedat,edgethn,'same');
edgedat = edgedat>0;
edgedat(braindat > 0.0) = 0; 

Vedge = spm_vol(funcmask);
Vedge.fname = fullfile(aroma_dir, 'edgedata.nii'); 
Vedge.descrip = 'edgedata';
Vedge = rmfield(Vedge, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
spm_write_vol(Vedge, edgedat);

icaparams_file = fullfile(ica_dir,'ica_aroma_ica_parameter_info.mat');

load(icaparams_file);

%%%%%%%%%% Get the required variables from sesInfo structure %%%%%%%%%%
% Number of subjects
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

% Number of components
numComp = sesInfo.numComp;

dataType = sesInfo.dataType;

mask_ind = sesInfo.mask_ind;

% First scan
structFile = deblank(sesInfo.inputFiles(1).name(1, :));

%%%%%%%% End for getting the required vars from sesInfo %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Get component data %%%%%%%%%%%%%%%%%

% Get the ICA Output files
icaOutputFiles = sesInfo.icaOutputFiles;

[subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', icaOutputFiles, 'numOfSub', ...
        numOfSub, 'numOfSess', numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);

% component files
if ~exist(meanALL_ICAFile.name)
    compFiles = subjectICAFiles(1).ses(1).name;
else
    compFiles = meanALL_ICAFile.name;
end

compFiles = icatb_fullFile('directory', ica_dir, 'files', compFiles);

[compData, icaTimecourse, structuralImage, coords, HInfo, text_left_right, classComp] = icatb_loadICAData('structfile',structFile,'compfiles',compFiles,'comp_numbers',[1:numComp],...
                                                                                                          'converttoz','yes','threshold', 1,'returnvalue',3,'open_dialog','no');

% Reshape icasig to components by voxels
compData = permute(compData, [2 3 4 1]);

% Structural volume
dim = HInfo.DIM;
tdim = size(icaTimecourse);

% Reshape to 2d
compData = reshape(compData, [prod(dim(1:3)),numComp]);

%       """Fraction component outside GM or WM""

nobrainFract = zeros(numComp,1);
brainFract = zeros(numComp,1);
edgeFract = zeros(numComp,1);

for i = 1:numComp       
    
    nbmaskdat = nobrain;
    nbCompdat = compData(:,i);

    nbmaskdat = reshape(nbmaskdat,[numel(nbmaskdat),1]);
    nbCompdat = reshape(nbCompdat,[numel(nbmaskdat),1]);

    [comparison, ~, ~, ~] = icatb_correlateFunctions(nbmaskdat, nbCompdat);
    nobrainFract(i) = comparison;

    brainmaskdat = braindat;
    brainCompdat = compData(:,i);

    brainmaskdat = reshape(brainmaskdat,[numel(brainmaskdat),1]);
    brainCompdat = reshape(brainCompdat,[numel(brainmaskdat),1]);

    [comparison, ~, ~, ~] = icatb_correlateFunctions(brainmaskdat, brainCompdat);
    brainFract(i) = comparison;

    edgemaskdat = edgedat;
    edgeCompdat = compData(:,i);

    edgemaskdat = reshape(edgemaskdat,[numel(edgemaskdat),1]);
    edgeCompdat = reshape(edgeCompdat,[numel(edgemaskdat),1]);

    [comparison, ~, ~, ~] = icatb_correlateFunctions(edgemaskdat, edgeCompdat);
    edgeFract(i) = comparison;
end

FT = abs(fft(icaTimecourse, [], 1));% FT of each IC along the firt dimension (time) 
FT = FT(1:(length(FT)/2) +1,:); % keep postivie frequencies (Hermitian symmetric) +1 to get Nyquist frequency bin as well 

%    """High frequency content"""
       
%   Determine sample frequency
    Fs = 1/t_r;
    
%   Determine Nyquist-frequency
    Ny = Fs/2;
    
%   Determine which frequencies are associated with every row in the melodic_FTmix file  (assuming the rows range from 0Hz to Nyquist)
    f = Ny * (1 : size(FT, 1)) / size(FT, 1);
    
%   Only include frequencies higher than 0.01Hz
    fincl = find(f > 0.01); %get indices
    FT = FT(fincl, :);
    f = f(fincl);
    
%     Set frequency range to [0-1]
    f_norm = (f - 0.01) / (Ny - 0.01);
    
%     For every IC; get the cumulative sum as a fraction of the total sum
    fcumsum_fract = cumsum(FT,1) ./ sum(FT,1);
    
%     Determine the index of the frequency with the fractional cumulative sum closest to 0.5
    [~, idx_cutoff] = min(abs(fcumsum_fract - 0.5));
    
%     Now get the fractions associated with those indices, these are the final feature scores
    HFC = f_norm(idx_cutoff)';
    
%         """Maximum robust correlation with confounds"""  
    
%     Read confounds
    conf_model = load(ppparams.rp_file); 

%     Determine the maximum correlation between confounds and IC time-series
    [nmixrows, nmixcols] = size(icaTimecourse);
    [nconfrows, nconfcols] = size(conf_model);
    
%     Resample the mix and conf_model matrices to have the same number of columns
    nmix = repmat(icaTimecourse, 1, compute_lcm(nmixcols,nconfcols)/nmixcols); 
    nconf_model = repmat(conf_model, 1, compute_lcm(nmixcols,nconfcols)/nconfcols);
    
    corr_mat = corr(nmix, nconf_model);
    max_correls = max(corr_mat, [], 2);
    
    max_correls = double(reshape(max_correls,[nmixcols,int64(size(max_correls, 1)/nmixcols)]));
    maxRPcorr = max_correls(:,1);

%      """ This function classifies a set of components into motion and non-motion components based on three features; 
%     maximum RP correlation, high-frequency content and non-brain-fraction"""
    
%   Define criteria needed for classification (thresholds and hyperplane-parameters)
    thr_corr = 0.50;
    thr_HFC = f_norm(min(find((f-0.12)>0)));

%   Classify the ICs

    noiseICs = find((nobrainFract > brainFract) | (edgeFract > brainFract) | (HFC > thr_HFC) | (maxRPcorr > thr_corr));

    noiseICdata = icaTimecourse(:,noiseICs);
    
%   Save results

    varnames = {'NoBrain_fraction','Edge-fraction','Brain_fraction','high_frequency_content','max_correlations'};
    T = table(nobrainFract,edgeFract,brainFract,HFC,maxRPcorr,'VariableNames',varnames);

    aroma_file = fullfile(aroma_dir,'AROMA_desission.csv');

    writetable(T,aroma_file,'WriteRowNames',false);
           
    [funcpath,funcfname,~] = fileparts(ppparams.funcfile{1});

    fg = spm_figure('FindWin','Graphics');

    for icomp=1:numComp
        spm_figure('Clear','Graphics');
        plot([0:t_r:t_r*(tdim(1)-1)],icaTimecourse(:,icomp))
        
        saveas(fg,fullfile(aroma_dir,['comp-' num2str(icomp,'%03d') '_time.png']));

        spm_figure('Clear','Graphics');
        plot(f,FT(:,icomp)./max(FT(:,icomp),[],'all'),f,fcumsum_fract(:,icomp))
        xline(thr_HFC*Ny)
        
        saveas(fg,fullfile(aroma_dir,['comp-' num2str(icomp,'%03d') '_frequency.png']));
    end
    
    % Path to noise ICs (might want to change path)
    noise_ICs_dir = fullfile(aroma_dir, strcat('ICA-AROMA_ICs_noise_',funcfname,'.txt'));
    
    if ~isempty(noiseICs) 
        if length(size(noiseICs)) > 0
            dlmwrite(noise_ICs_dir, noiseICs(:), 'precision', "%i"); %write matrix to text file with values separated by a ',' as a delimiter
        else
            dlmwrite(noise_ICs_dir, int64(noiseICs), 'precision', "%i");
        end
    else
        f = open(confounds_ICs,'w');
        fclose(f);
    end

end
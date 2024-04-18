function acid_defaults

% Created at 2020-10-19

global acid_def;

acid_def.local_defaults = {fullfile(fileparts(mfilename('fullpath')),'local','acid_local_defaults.m')};


%% GENERAL

% Tolerance for rounding of bvalues (bval_rounded = round(bvals/rounding_tolerance)*rounding_tolerance;)
acid_def.bval.rounding_tolerance = 100;


%% ---RESAMPLE---

% Interpolation order as defined in spm_reslice. A splin-interpolation is used by default.
acid_def.resample.interpol_reslice = -7;

%% ---Reliability masking

% Lowest tested model-fit error.
acid_def.relmask.relmask_thr_errmin = 1;

% Highest tested model-fit error.
acid_def.relmask.relmask_thr_errmax = 3;

% Number of steps.
acid_def.relmask.relmask_thr_nsteps = 100;

%% ---BRAIN MASK---

% This option reduces holes in the brain mask by smoothing the mask or the tissue probability maps. If it is set to zero no smoothing is used.
acid_def.brainmask.smk   = [3 3 3];

% Brain coverage (0: brain covers full field-of-view).
acid_def.brainmask.perc1 = 0.8;

%% ---DTI/DKI---

% Write fitted diffusion-weighted images (DWI). This option provides DWIs with reduced noise, but keep in mind that the modeled DWIs are limited by the DTI model used.
acid_def.diffusion.dummy_DTIsm = 0;

% If the vendor uses another coordinate system than the coordinate system, in which your b-vectors were defined, you need to reorient them. Provide a 3 x 3  - matrix to reorient b-vectors.
% The b-vectors are rotated to be compliant with FSL.
acid_def.diffusion.rmatrix = [1  0  0;
                              0 -1  0; 
                              0  0  1];


acid_def.diffusion.dummy_plot_error = 0;

%%% OUTPUTS

% Write Freiburg Tractography format. The path of the Freiburg Tractography
% tools have to be included - otherwise it will not work. 
% IMPORTANT: The rmatrix (see above) has to be changed to [1  0  0;
%                                                          0  1  0; 
%                                                          0  0  1];
acid_def.diffusion.dummy_write_TFreiburg = 0;

% dummy for writing out 4d model-fit error. This option provides DWIs with reduced noise, but keep in mind that the modeled DWIs are limited by the DTI model used.
acid_def.diffusion.dummy_write_error = 0;

% dummy for writing out 4d model-fit error per shell (only for NLLS). This option provides DWIs with reduced noise, but keep in mind that the modeled DWIs are limited by the DTI model used.
acid_def.diffusion.dummy_error_per_shell = 0;

% Option for OLS to write out non log RMSE. 0 == RMSE-LOG , 1 == RMSE.
acid_def.diffusion.dummy_log = 0;

% dummy for writing out 4d fitted signal. This option provides DWIs with reduced noise, but keep in mind that the modeled DWIs are limited by the DTI model used.
acid_def.diffusion.dummy_write_fitted = 0;

% dummy for writing out all eigenvalues and eigenvectors
acid_def.diffusion.dummy_eigensys = 1;

% dummy for writing out tensors
acid_def.diffusion.dummy_write_tensors = 0;

% WLS

% Option to save the weights that were generated during robust fitting on a voxel-by-voxel basis as a 4d volume.
acid_def.diffusion.dummy_weights_WLS = 0;

% ROBUST FITTING

% Option to save the weights that were generated during robust fitting on a voxel-by-voxel basis as a 4d volume.
acid_def.diffusion.dummy_weights_ROBUST = 0;

% NLLS

% Option to save the uncertainty map.
acid_def.diffusion.dummy_calculate_uncertainty = 0;

%%%

% dummy for masking output
% if the value is set to 0 and a mask is still defined, the output won't be masked, BUT:
% - the diagnostic plots of model-fit error are computed only wihtin the mask
% - the weights are computed only within the mask (applicable to the slice-wise weight in robust fitting)
acid_def.diffusion.dummy_mask = 1;

% dummy for applying strict conditions on the DTI maps. If set to YES, we
% require positive definitess of the diffusion tensor
acid_def.diffusion.dummy_cond_strict = 0;

% This is a the value of a Tikhonov regularization constant and should not exceed 1e-4 otherwise it can bias the tensor estimates. It is also related to the number of iteration taking place.

acid_def.diffusion.thr_DTvar = 1e-4;

% Condition number theshold

acid_def.diffusion.thr_cond = 1e+5;

% This option reduces holes in the brain mask by smoothing the mask or the tissue probability maps. If it is set to zero no smoothing is used.
acid_def.diffusion.smk     = [3 3 3];

% Confidence interval for robust fitting.
acid_def.diffusion.cval    = 0.3;

% Threshold for minimal diffusivity.
acid_def.diffusion.dthr    = 1e-7;

% This variable determines the smoothing that is applied on the residuals; smaller maxk results in more smoothing.
acid_def.diffusion.kmax    = 16;

% Standard deviation of logarithm of the signal outside the brain. This measure is used if the noise cannot be estimated from outside the brain (see brain mask option).
acid_def.diffusion.sigma0  = 4;

% Determines whether to write out all eigenvectors of the diffusion tensor, i.e. 1st, 2nd, and 3rd eigenvector and eigenvalues.The filename will be extended by the number of the eigenvector componenten and the eigenvector number, i.e. "filename-ij.img" with "i" being the eigenvector componenten and "j" being the eigenvector number.;
acid_def.diffusion.dummyDT = 0;

% Maximum number of iteration in the IRLS approach underlying robust tensor fitting. Note that the iteration terminates earlier if the rms difference between the previous and actual tensor falls below the Tikhonov regularization factor specified above.
acid_def.diffusion.Niter = 10;

%%% Robust Fitting

% Standard deviation of logarithm of the signal outside the brain. This measure is used if the noise cannot be estimated from outside the brain (see brain mask option).
acid_def.diffusion.robust_sigma = '';

% A1 scaling factor for w1. Scaling factor for regional weight. The value is set to 0.3 (smaller C means less outliers will be rejected, see Zwiers et al., Neuroimage 2010).
acid_def.diffusion.A1 = 0.3;

% A2 scaling factor for w2. Scaling factor for slice-wise weight. The value is set to 0.3 (smaller C means less outliers will be rejected, see Zwiers et al., Neuroimage 2010).
acid_def.diffusion.A2 = 0.1;


% Option to display weights during robust fitting. Might be instructive to look at it once in a while.
acid_def.diffusion.dummy_plot = 0;

%%% NLLS

% Choose the number of fit iterations.
acid_def.diffusion.max_iter_GN = 75;

% Number of coils used in your measurement, the noise estimate (sigma) must have been done with the same number of coils L.
acid_def.diffusion.n_coils = 1;








%% ---ECMO---

% Choose whether you want to see the estimated EC and motion parameter for each image.
acid_def.ecmo.plot = 0;

% Choose whether you want to reslice and output the registered images. (default: ON)
acid_def.ecmo.write = 1;

% If "true", the iteration for the transformation parameters (TP) is initialized by the TPs of the b0 volumes.'
%     'It is the recommended option if the b0 and higher b-shells are acquired interspersed.' Can be "true" or "false".
acid_def.ecmo.initialization = true;

% Separation for resampling
acid_def.ecmo.separation = [4 2];

% Interpolation order as defined in spm_reslice. A splin-interpolation is used by default.
acid_def.ecmo.interpol_reslice = -7;

% Choose the parameters, which you want to correct. You can choose between 12 affine parameters. The 4 eddy current parameters are displayed in Figure 1 (see ECMOCO_Documentation). We propose three sets of parameters for different purposes (see below), but you can select the parameters freely.
%acid_def.ecmo.freeze = [ones(1,6), 0, 1, 0, 1, 1, 0];

% Averaging in z-direction (only for slice-wise registration). Number of slices for smoothing in z-direction can be set. 1 == deactivated'
acid_def.ecmo.zsmooth  = 1;


%% ---HYSCO---

% Choose regularization parameter alpha [x,y,z] that weights the ''diffusion'' regularizer. For larger values of alpha, the computed solution will in general be smoother, but the image distance between the corrected blip-up and blip-down image will be larger.
acid_def.hysco.alpha     = [50 50 50];

% Choose regularization parameter beta that weights the ''Jacobian'' regularizer that guarantees the invertibility of the transformations. By design the value of this regularization functional grows to infinity when the Jacobian determinant of either of the geometrical transformations approaches zero. Thus, for larger/smaller values of beta, the range of the Jacobian determinants, which translates to the maximum compression/expansion of volume, becomes smaller/larger. However, for any beta that is greater than zero both transformations will be invertible.
acid_def.hysco.beta      = 10;

% Apply the field inhomogeneity estimated from the reference images to the
% other image volumes'' (see 3 and 4). If set to ''no'' and if the same number of diffusion-weighted images is provided for blip-up and blip-down the susceptibility correction is done for each image separately (This might be useful to correct for the distortions induced by nonlinear eddy current fields). To this end, the field-inhomogeneity estimated from the non-diffusion weighted images is used as a starting guess for minimization of the distance between the respective diffusion-weighted image pairs. Optimization is only carried out on the finest discretization level to save computation time.
acid_def.hysco.dummy_ecc = 0;

%{
Choose parameter theta that balances between the data fit and the smoothness of the cubic B-spline approximation of the image data. For theta equal to zero a standard cubic B-spline interpolation is used. For positive theta data is only approximated, but the image representation is smoother. Thus, theta can be used to adjust for noise level of the data. Note that this scheme is only used for the optimization and in particular that the corrected image data is obtained by resampling using a standard tri-linear interpolation.
For details see Section 3.4 in:'
Modersitzki, J. FAIR: Flexible Algorithms for Image Registration. Society for Industrial and Applied Mathematics (SIAM); 2009.
%}
acid_def.hysco.theta = 1e-1;

% NOT RECOMMENDED! (1 == No ; 0 == Yes) Apply the field inhomogeneity estimated from the reference images to the ''other image volumes'' (see 3 and 4). If set to ''no'' and if the same number of diffusion-weighted images is provided for blip-up and blip-down the susceptibility correction is done for each image separately (This might be useful to correct for the distortions induced by nonlinear eddy current fields). To this end, the field-inhomogeneity estimated from the non-diffusion weighted images is used as a starting guess for minimization of the distance between the respective diffusion-weighted image pairs. Optimization is only carried out on the finest discretization level to save computation time.
acid_def.hysco.dummy_ecc = 0;

% Interpolation order as defined in spm_reslice. A splin-interpolation is used by default.
acid_def.hysco.resample = -7;

%{
Choose data dimensions that are coarsened in multilevel optimization.

Background: To accelerate computations and improve robustness against local minima, we use a multilevel
strategy for optimization. One part of this is to represent the image data by coarsened versions using averaging.
This setting affects which dimensions are coarsened.

Default: [1,1,1] coarsens data along all coordinate directions

Example: choosing [1,1,0] avoids coarsening along the third direction (commonly, the slice selection direction).
This choice is motivated by the fact that EPI distortions are typically restricted within slices.

Note: This flag does not affect the final resolution of the estimated inhomogeneity, which will be
related to the number of voxels in the original image data.
%}
acid_def.hysco.restrictdim = [1 1 1];




end

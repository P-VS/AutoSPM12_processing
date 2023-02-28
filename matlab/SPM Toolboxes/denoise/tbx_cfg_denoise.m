function denoise = tbx_cfg_denoise

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','denoise')); end

% Regression
%---------------------------------------
func_file              = cfg_files;
func_file.tag          = 'func_file';
func_file.name         = 'Functional data';
func_file.filter       = 'image';
func_file.ufilter      = '.*';
func_file.num          = [1 Inf];
func_file.help         = {'Select the fMRI data'};

rp_file              = cfg_files;
rp_file .tag          = 'rp_file';
rp_file .name         = 'Motion parameter file';
rp_file .filter       = 'mat';
rp_file .ufilter      = '.*';
rp_file .num          = [1 1];
rp_file .help         = {'Select the motion parameter file (rp_... file)'};

expand_regressors         = cfg_menu;
expand_regressors.tag     = 'expand_regressors';
expand_regressors.name    = 'Regressors expansion';
expand_regressors.help    = {['Expand the motion regressors? ' ...
                              'No (6 regressors) / Derivatives (12 regressors) / Squared (24 regressors)']};
expand_regressors.values  = {1 2 3};
expand_regressors.val     = {1};
expand_regressors.labels  = {'No' 'Derivatives' 'Squared'};

% aCompCor
%---------------------------------------
csf_map              = cfg_files;
csf_map.tag          = 'csf_map';
csf_map.name         = 'CSF map';
csf_map.filter       = 'image';
csf_map.ufilter      = '.*';
csf_map.num          = [1 1];
csf_map.help         = {'Select the CSF map (c3...)'};

% Main structure
%---------------------------------------
acompcor           = cfg_exbranch;
acompcor.tag       = 'acompcor';
acompcor.name      = 'aCompCor';
acompcor.help      = {['This function denoises fMRI data using aCompCor. '...
                       'Those components that explain 50% of the signal variation in CSF will be removed.']};
acompcor.val       = {func_file csf_map rp_file expand_regressors csf_map};
acompcor.prog      = @(job)vout_denoise('acompcor',job);
acompcor.vout      = @(job)vout_denoise('vout',job);

regression           = cfg_exbranch;
regression.tag       = 'regression';
regression.name      = 'Regression';
regression.help      = {'This function denoises fMRI data using motion regression.'};
regression.val       = {func_file rp_file expand_regressors};
regression.prog      = @(job)vout_denoise('regression',job);
regression.vout      = @(job)vout_denoise('vout',job);

denoise         = cfg_choice;
denoise.tag     = 'fmri_denoise';
denoise.name    = 'fMRI denoising';
denoise.help    = {'This toolbox is meant to denoise fMRI data'};
denoise.values  = {regression acompcor};

function out = vout_denoise(cmd,job)

switch lower(cmd)
    case 'regression'
        [out.dnfiles out.norms]=denoise_regression(job);
    case 'acompcor'
        [out.dnfiles out.norms]=denoise_acompcor(job);
    case 'vout'
        out(1)           =cfg_dep;
        out(1).sname     =sprintf('Denoised fMRI files');
        out(1).src_output=substruct('.','dnfiles');
        out(1).tgt_spec  =cfg_findspec({{'filter','image','strtype','e'}});
        
        out(2)           =cfg_dep;
        out(2).sname     =sprintf('Denoised rp_file');
        out(2).src_output=substruct('.','norms');
        out(2).tgt_spec  =cfg_findspec({{'filter','rp_*','strtype','e'}});
end

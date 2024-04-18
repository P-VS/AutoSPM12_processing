function acid_config = tbx_main_acid_config
% Configuration file for the "histological MRI" (acid) toolbox
% Previously named "Voxel-Based Quantification" (VBQ)
% -> Dealing with local defaults ("Configure Toolbox" branch)
%_______________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging


%--------------------------------------------------------------------------
% Customised defaults
%--------------------------------------------------------------------------
customised         = cfg_files;
customised.tag     = 'customised';
customised.name    = 'Customised';
customised.help    = {['Select the [acid_local_defaults_*.m] file containing ' ...
    'the specific defaults to process your data. Note that all other defaults ' ...
    'values will be reinitialised to their standard values.']};
customised.filter  = 'm';
customised.dir     = fullfile(fileparts(mfilename('fullpath')),'local');
customised.ufilter = '^acid_.*\.m$';
customised.num     = [1 1];
%customised.def     = @(val)acid_get_defaults('local_defaults', val{:});

% ---------------------------------------------------------------------
% Use standard defaults (no customization)
% ---------------------------------------------------------------------
standard           = cfg_entry;
standard.tag       = 'standard';
standard.name      = 'Standard';
standard.help      = {''};
standard.strtype = 's';
standard.num     = [1 Inf];
standard.val     = {'yes'};

%--------------------------------------------------------------------------
% Set acid defaults parameters
%--------------------------------------------------------------------------
acid_setdef         = cfg_choice;
acid_setdef.tag     = 'acid_setdef';
acid_setdef.name    = 'Defaults parameters';
acid_setdef.help    = {['You can either stick with standard defaults parameters ' ...
    'from [acid_defaults.m] or select your own customised defaults file.']};
acid_setdef.values  = {standard customised};
acid_setdef.val     = {standard};


% ---------------------------------------------------------------------
% Configure the acid Toolbox - load local, user-defined defaults file and
% overwrite standard defaults 
% ---------------------------------------------------------------------
acid_config         = cfg_exbranch;
acid_config.tag     = 'acid_config';
acid_config.name    = 'Configure toolbox';
acid_config.val     = { acid_setdef };
acid_config.help    = {['Customised default parameters can be set here by selecting ' ...
    'a customised [acid_local_defaults_*.m] file. Type [help acid_local_defaults] for ' ...
    'more details.']};
acid_config.prog    = @acid_run_config;

end
%----------------------------------------------------------------------

% =========================================================================
% (VOUT &) RUN SUBFUNCTION(S)
% =========================================================================
function out = acid_run_config(job)
%==========================================================================
% PURPOSE
% Load standard defaults and overwrite thme by customised values if
% provided via local defaults file. 
%==========================================================================

% reinitialise standard defaults
acid_defaults;
% overwrite with customised defaults
if isfield(job.acid_setdef,'customised')
    deffnam = job.acid_setdef.customised;
    acid_get_defaults('local_defaults',deffnam);
    spm('Run',deffnam);
end
out = acid_get_defaults;

end
%_______________________________________________________________________

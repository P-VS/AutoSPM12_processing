function do_ica(mask, t_r, ppparams)
% DO_ICA Creates the parameter file necessary to run ICA with the GIFT toolbox
% and performs ICA

% mask = func mask
% funcdata = func data
% t_r = repetition time (read from json file)

if numel(ppparams.echoes)>1 && ~ppparams.mecombined
    n_sessions = numel(ppparams.echoes);
else
    n_sessions = 1;
end

%% Estimate number of components (MDL - FWHM)
dim_est_opts.method = 2; %MDL
dim_est_opts.fwhm = [5,5,5]; %FWHM
comp_est = 0;

for ie=ppparams.echoes
    [tcomp_est, mdl, aic, kic] = icatb_estimate_dimension(fullfile(ppparams.ppfuncdir,ppparams.func(ie).sfuncfile), mask, 'double', dim_est_opts);  

    comp_est = comp_est + tcomp_est;
end

comp_est = round(comp_est / n_sessions);

%% Adapt template parameter to specific subject data

icatb_defaults;

% Load template parameters .mat file
icaparms_file = ica_aroma_make_gift_parameters(t_r,n_sessions,mask,comp_est,ppparams);

%% Set up ICA
icaparam_file = icatb_setup_analysis(icaparms_file);

load(icaparam_file);

%% Run Analysis (All steps)
sesInfo = icatb_runAnalysis(sesInfo, 1);

end



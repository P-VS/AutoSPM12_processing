function relmask_thr = tbx_cfg_acid_relmask_thr

% model-fit error maps
ERROR         = cfg_files;
ERROR.tag     = 'ERROR';
ERROR.name    = 'Model-fit error maps';
ERROR.help    = {'Select the model-fit error maps.'};
ERROR.filter  = 'image';
ERROR.ufilter = '.*';
ERROR.num     = [0 Inf];

% FA maps
FA         = cfg_files;
FA.tag     = 'FA';
FA.name    = 'FA maps';
FA.help    = {'Select the FA maps.'};
FA.filter  = 'image';
FA.ufilter = '.*';
FA.num     = [0 Inf];

% subject-specific masks
maskSUB         = cfg_files;
maskSUB.tag     = 'maskSUB';
maskSUB.name    = 'Subject-specific masks';
maskSUB.help    = {'Select subject-specific masks.'};
maskSUB.filter  = 'image';
maskSUB.ufilter = '.*';
maskSUB.num     = [0 Inf];

% global mask
maskGLO         = cfg_files;
maskGLO.tag     = 'maskGLO';
maskGLO.name    = 'Global mask';
maskGLO.help    = {'Select a global mask.'};
maskGLO.filter  = 'image';
maskGLO.ufilter = '.*';
maskGLO.num     = [0 Inf];
%{
% start
relmask_thr_errmin         = cfg_entry;
relmask_thr_errmin.tag     = 'relmask_thr_errmin';
relmask_thr_errmin.name    = 'Lowest tested model-fit error.';
relmask_thr_errmin.help    = {'Lowest tested model-fit error.'};
relmask_thr_errmin.strtype = 'r';
relmask_thr_errmin.num     = [1 1];
relmask_thr_errmin.val     = {1};

% end
relmask_thr_errmax         = cfg_entry;
relmask_thr_errmax.tag     = 'relmask_thr_errmax';
relmask_thr_errmax.name    = 'Highest tested model-fit error.';
relmask_thr_errmax.help    = {'Highest tested model-fit error.'};
relmask_thr_errmax.strtype = 'r';
relmask_thr_errmax.num     = [1 1];
relmask_thr_errmax.val     = {3};

% number of steps
relmask_thr_nsteps          = cfg_entry;
relmask_thr_nsteps.tag      = 'relmask_thr_nsteps';
relmask_thr_nsteps.name    = 'Number of steps.';
relmask_thr_nsteps.help    = {'Number of steps.'};
relmask_thr_nsteps.strtype = 'r';
relmask_thr_nsteps.num     = [1 1];
relmask_thr_nsteps.val     = {100};
%}
% display figure
relmask_thr_plot           = cfg_menu;
relmask_thr_plot.tag       = 'relmask_thr_plot';
relmask_thr_plot.name      = 'Choose display option';
relmask_thr_plot.help      = {'Choose whether you want to see the plot.'
                       'default: ON'
                         };
relmask_thr_plot.labels    = {
                                'ON'
                                'OFF'
                                    }';
relmask_thr_plot.values    = {1 0};
relmask_thr_plot.val       = {1};

% exbranch
relmask_thr      = cfg_exbranch;
relmask_thr.tag  = 'relmask_thr';
relmask_thr.name = 'Determine threshold';
relmask_thr.help = {'Post-processing'};
relmask_thr.val  = {ERROR, FA, maskSUB, maskGLO, relmask_thr_plot};
relmask_thr.prog = @local_REL_determine_thr;
relmask_thr.vout = @vout_REL_determine_thr;

end

function out = local_REL_determine_thr(job)

    relmask_thr_errmin = acid_get_defaults('relmask.relmask_thr_errmin');
    relmask_thr_errmax = acid_get_defaults('relmask.relmask_thr_errmax');
    relmask_thr_nsteps = acid_get_defaults('relmask.relmask_thr_nsteps');

    thr = acid_relmask_thr(char(job.ERROR), char(job.FA), char(job.maskSUB), ...
        char(job.maskGLO), relmask_thr_errmin, relmask_thr_errmax,...
        relmask_thr_nsteps, job.relmask_thr_plot);
    
    out.thr = thr;
end

function dep = vout_REL_determine_thr(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Threshold for reliability masking';
    dep(1).src_output = substruct('.','thr');
    dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});
end
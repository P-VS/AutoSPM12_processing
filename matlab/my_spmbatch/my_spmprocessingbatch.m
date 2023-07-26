function matlabbatch = my_spmprocessingbatch(sub,ses,task,params)

substring = ['sub-' num2str(sub,'%02d')];
subpath = fullfile(params.datpath,substring,['ses-' num2str(ses,'%03d')]);

subfmridir = fullfile(subpath,'func');

if params.multi_echo && params.use_echoes_as_sessions
    nechoes = numel(params.echoes);
else
    nechoes = 1;
end

step = 0;
for k=1:numel(task)
    eventsdat{k} = fullfile(subfmridir,[substring '_task-' task{k} '_events.tsv']);

    if params.multi_echo
        jsonfile{k} = fullfile(subfmridir,[substring '_task-' task{k} '_bold_e1.json']);
    else
        jsonfile{k} = fullfile(subfmridir,[substring '_task-' task{k} '_bold.json']);
    end

    for ne=1:nechoes
        if params.multi_echo && params.use_echoes_as_sessions
            endfix = [params.fmri_endfix '_e' num2str(params.echoes(ne))];
        else
            endfix = params.fmri_endfix;
        end

        funcfile = fullfile(subpath,params.preprocfmridir,[params.fmri_prefix substring '_task-' task{k} '_' endfix '.nii']);
        Vfunc = spm_vol(funcfile);
    
        for i=1:numel(Vfunc)
            ppfmridat{k}.sess{ne}.func{i,1} = [Vfunc(i).fname ',' num2str(i)];
        end
    end

    if params.add_regressors
        confound_file{k} = fullfile(subpath,params.preprocfmridir,[params.confounds_prefix substring '_task-' task{k} '_bold.txt']);
    end
end

%fMRI model specification
resultmap = fullfile(subpath,['SPMMAT-' task{k} '_' params.analysisname]);
if exist(resultmap,'dir'); rmdir(resultmap,'s'); end
mkdir(resultmap)

jsondat = fileread(jsonfile{1});
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;
SliceTiming = jsondat.SliceTiming;
nsl= ceil(numel(SliceTiming)/numel(find(SliceTiming==SliceTiming(1))));

step=step+1;
matlabbatch{step}.spm.stats.fmri_spec.dir = {resultmap};
matlabbatch{step}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{step}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{step}.spm.stats.fmri_spec.timing.fmri_t = nsl;
matlabbatch{step}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

for k=1:numel(task)
    %correct events file for dummy scans if needed
    dummys = floor(params.dummytime/tr);
    edat{k} = tdfread(eventsdat{k},'\t');
    edat{k}.onset = edat{k}.onset-dummys*tr;

    for it=1:numel(edat{k}.trial_type(:,1))
        ntrial_type(it,1) = convertCharsToStrings(edat{k}.trial_type(it,:));
    end

    edat{k}.trial_type = ntrial_type;
    
    [~,edatorder] = sort(edat{k}.trial_type);
    edat{k}.onset = edat{k}.onset(edatorder);
    edat{k}.duration = edat{k}.duration(edatorder);
    edat{k}.trial_type = edat{k}.trial_type(edatorder,:);

    numc=0;
    for trial=1:numel(edat{k}.onset)
        if isstring(edat{k}.trial_type(trial,:))
            trial_type = edat{k}.trial_type(trial,:);
        else
            trial_type = num2str(edat{k}.trial_type(trial,:));
        end

        if numc>0
            numc = numel(edat{k}.conditions)+1;
            for nc=1:numel(edat{k}.conditions)
                if strcmp(edat{k}.conditions{nc}.name,strtrim(trial_type)); numc=nc; end
            end
            if numc<numel(edat{k}.conditions)+1
                edat{k}.conditions{numc}.onsets = [edat{k}.conditions{numc}.onsets edat{k}.onset(trial)];
                edat{k}.conditions{numc}.durations = [edat{k}.conditions{numc}.durations edat{k}.duration(trial)];
            else
                edat{k}.conditions{numc}.name = strtrim(trial_type);
                edat{k}.conditions{numc}.onsets = [edat{k}.onset(trial)];
                edat{k}.conditions{numc}.durations = [edat{k}.duration(trial)];
            end
        else
            edat{k}.conditions{1}.name = strtrim(trial_type);
            edat{k}.conditions{1}.onsets = edat{k}.onset(trial);
            edat{k}.conditions{1}.durations = edat{k}.duration(trial);
            numc=1;
        end
    end

    for ne=1:nechoes
        nsess = (k-1)*nechoes+ne;

        matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).scans = ppfmridat{k}.sess{ne}.func(:,1);

        for nc=1:numel(edat{k}.conditions)
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).cond(nc).name = char(edat{k}.conditions{nc}.name);
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).cond(nc).onset = edat{k}.conditions{nc}.onsets;
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).cond(nc).duration = edat{k}.conditions{nc}.durations;
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).cond(nc).tmod = 0;
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).cond(nc).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).cond(nc).orth = 1;
        end

        matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).multi = {''};
        matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).regress = struct('name', {}, 'val', {});
    
        if params.add_regressors
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).multi_reg = {confound_file{k}};
        else
            matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).multi_reg = {''};
        end
        matlabbatch{step}.spm.stats.fmri_spec.sess(nsess).hpf = params.hpf;
    end
end

matlabbatch{step}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{step}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{step}.spm.stats.fmri_spec.volt = 1;
matlabbatch{step}.spm.stats.fmri_spec.global = 'None';

matlabbatch{step}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{step}.spm.stats.fmri_spec.mask = {''};

matlabbatch{step}.spm.stats.fmri_spec.cvi = params.model_serial_correlations;

step=step+1;
%Model estimation
matlabbatch{step}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{step-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{step}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{step}.spm.stats.fmri_est.method.Classical = 1;

step=step+1;
%Contrast Manager
matlabbatch{step}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{step-1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
for ic=1:numel(params.contrast)
    contrastname='';
    weights = [];
    for k=1:numel(task)
        subweights = zeros(1,numel(edat{k}.conditions));
        if params.add_regressors
            rpdat = load(confound_file{k});
            subweights=[subweights zeros(1,numel(rpdat(1,:)))]; 
        end
        for icn=1:numel(params.contrast(ic).conditions)
            if k==1
                if params.contrast(ic).vector(icn)>0; contrastname = [contrastname ' + ' params.contrast(ic).conditions{icn}]; end
                if params.contrast(ic).vector(icn)<0; contrastname = [contrastname ' - ' params.contrast(ic).conditions{icn}]; end
            end
            
            indx=0;
            for icn2=1:numel(edat{k}.conditions)
                if strcmp(lower(params.contrast(ic).conditions{icn}),lower(edat{k}.conditions{icn2}.name)); indx=icn2; end
            end

            if indx>0; subweights(indx)=params.contrast(ic).vector(icn); end
        end
        weights=[weights subweights];
    end

    matlabbatch{step}.spm.stats.con.consess{ic}.tcon.name = contrastname;
    matlabbatch{step}.spm.stats.con.consess{ic}.tcon.weights = weights;

    if params.multi_echo && params.use_echoes_as_sessions
        matlabbatch{step}.spm.stats.con.consess{ic}.tcon.sessrep = 'replsc';
    else
        matlabbatch{step}.spm.stats.con.consess{ic}.tcon.sessrep = 'none';
    end
end
matlabbatch{step}.spm.stats.con.delete = 0;

end
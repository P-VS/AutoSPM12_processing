function sigma = tbx_cfg_acid_noise_estimate(bvals)

%% sigma estimation
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'Select the images in which the noise is to be estimated. You can load in either a set of 3d nifti or a single 4d nifti file.'};
images.filter  = 'image';
images.ufilter = '.*';
images.num     = [0 Inf];

mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Mask';
mask.help    = {'Enter a binary mask that defines a region in which the noise is to be estimated.'};
mask.filter  = 'image';
mask.ufilter = '.*';
mask.num     = [1 1];

ncoils         = cfg_entry;
ncoils.tag     = 'ncoils';
ncoils.name    = 'ncoils';
ncoils.help    = {'Effective number of receiver coils.'};
ncoils.strtype = 'e';
ncoils.num     = [1 1];
ncoils.val     = {1e0};

ncoils_choice        = cfg_choice;
ncoils_choice.tag    = 'ncoils_choice';
ncoils_choice.name   = 'Standard method';
ncoils_choice.help   = {''};
ncoils_choice.values = {ncoils};
ncoils_choice.val    = {ncoils};

% bvals.name = 'Repeated measures method';



% dummy for writing out W maps
dummy_shell        = cfg_menu;
dummy_shell.tag    = 'dummy_shell';
dummy_shell.name   = 'Shell for repeated measurement';
dummy_shell.help   = {'This option allows you to select if either the b0 images or the images of the highest shell will be used for repetitions.'};
dummy_shell.labels = {
    'b0 images'
    'Highest shell'
    }';
dummy_shell.values = {0 1};
dummy_shell.val    = {0};

repeated      = cfg_branch;
repeated.tag  = 'repeated';
repeated.name = 'Repeated measures method';
repeated.help = {''};
repeated.val  = {bvals, dummy_shell};

noise_dummy_type        = cfg_choice;
noise_dummy_type.tag    = 'noise_dummy_type';
noise_dummy_type.name   = 'Noises estimation method';
noise_dummy_type.help   = {' '};
noise_dummy_type.values = {repeated, ncoils_choice};
noise_dummy_type.val    = {repeated};



sigma      = cfg_exbranch;
sigma.tag  = 'sigma';
sigma.name = 'Noise estimation';
sigma.val  = {images mask noise_dummy_type};
sigma.help = {
    'Select the images you want to calculate sigma from from as well as the noise mask.'
    'Changeable default: '
    '- ncoils: Number of receiver coils. Effective number can be less!'
    };
sigma.prog = @local_sigma_est;
sigma.vout = @vout_sigma_est;

end

function out = local_sigma_est(job)

% read in bvals and bvecs
if isfield(job.noise_dummy_type,'repeated')
    if isfield(job.noise_dummy_type.repeated,'bvals') && isfield(job.noise_dummy_type.repeated.bvals,'bvals_file')

        job.files1 = job.noise_dummy_type.repeated.bvals.bvals_file{1,1};
        [~,~,e] = fileparts(job.files1);

        if ~isempty(job.files1)
            switch lower(e)
                case '.mat'
                    if size(cell2mat(struct2cell(load(job.files1))),1) == 1
                        job.bvals = cell2mat(struct2cell(load(job.files1)));
                    else
                        error('Unexpected file extension!')
                    end

                case '.nii'
                    error('Unexpected file extension!')
                otherwise
                    if size(dlmread(job.files1),1) == 1
                        job.bvals = dlmread(job.files1);
                    else
                        error('Unexpected file extension!')
                    end
            end
        else
            job.bvals = '';
        end

    elseif isfield(job.noise_dummy_type.repeated,'bvals') && isfield(job.noise_dummy_type.repeated.bvals,'bvals_exp')
        if ~isempty(job.noise_dummy_type.repeated.bvals.bvals_exp)
            job.bvals = job.noise_dummy_type.repeated.bvals.bvals_exp;
        else
            job.bvals = '';
        end
    end
end
if isfield(job.noise_dummy_type,'ncoils_choice')
    [sigma_est, sigma_est_avg] = acid_noise_estimate_std(char(job.images), char(job.mask), job.noise_dummy_type.ncoils_choice.ncoils);
else
    rounding_tolerance = acid_get_defaults('bval.rounding_tolerance');
    [sigma_est, ~, ~] = acid_noise_estimate_repmeas(char(job.images), char(job.mask), job.bvals, job.noise_dummy_type.repeated.dummy_shell,rounding_tolerance);
    sigma_est_avg = sigma_est;
end

out.sigma_est     = sigma_est;
out.sigma_est_avg = sigma_est_avg;

end

function dep = vout_sigma_est(~)
dep(1)            = cfg_dep;
dep(1).sname      = 'Noise estimates';
dep(1).src_output = substruct('.','sigma_est');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = 'Average noise estimate';
dep(2).src_output = substruct('.','sigma_est_avg');
dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
end
function startup = tbx_cfg_acid_startup(p_out, bvals_bvecs)
% B.Fricke 17/12/2021

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

config_path = fileparts(mfilename('fullpath'));
config_path = [config_path(1:end-9) 'config'];
customised.dir     = fullfile(config_path,'local');
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
acid_setdef.name    = 'Default parameters';
acid_setdef.help    = {['You can either stick with standard defaults parameters ' ...
    'from [acid_defaults.m] or select your own customised defaults file.']};
acid_setdef.values  = {standard customised};
acid_setdef.val     = {standard};


% images
images         = cfg_files;
images.tag     = 'images';
images.name    = 'dMRI dataset (3D NIfTIs)';
images.help    = {'Select the dMRI dataset (3D NIfTIs).'};
images.filter  = 'image';
images.ufilter = '.*';
images.num     = [0 Inf];

% file name
filename         = cfg_entry;
filename.tag     = 'filename';
filename.name    = 'Specify filename (without .nii)';
filename.help    = {'Output filename.'};
filename.strtype = 's';

% exbranch extract
extract       = cfg_branch;
extract.tag   = 'extract';
extract.name  = 'Extract bval/bvec files and 4D data';
extract.help  = {'Extracts and saves bval and bvec from the input NIfTI images. If the input is a set of 3D NIfTI images, it converts them into a single 4D NIfTI image (set filename A BIDS compliant filename can be set for the output 4D file.'};
extract.val   = {images, filename};


% images
images_existing         = cfg_files;
images_existing.tag     = 'images_existing';
images_existing.name    = 'dMRI dataset (3D/4D NIfTI)';
images_existing.help    = {'Select the dMRI dataset (3D/4D NIfTI).'};
images_existing.filter  = 'image';
images_existing.ufilter = '.*';
images_existing.num     = [0 Inf];


% exbranch existing
existing       = cfg_branch;
existing.tag   = 'existing';
existing.name  = 'Use existing bval/bvec files and 3D/4D data';
existing.help  = {'Extracts and saves bval and bvec from the input NIfTI images. If the input is a set of 3D NIfTI images, it converts them into a single 4D NIfTI image (set filename A BIDS compliant filename can be set for the output 4D file.'};
existing.val   = {images_existing, bvals_bvecs, filename};

% exbranch existing
existing_only_image       = cfg_branch;
existing_only_image.tag   = 'existing_only_image';
existing_only_image.name  = 'Use existing 3D data without bval/bvec files';
existing_only_image.help  = {'Extracts and saves bval and bvec from the input NIfTI images. If the input is a set of 3D NIfTI images, it converts them into a single 4D NIfTI image (set filename A BIDS compliant filename can be set for the output 4D file.'};
existing_only_image.val   = {images_existing, filename};

only         = cfg_entry;
only.tag     = 'only';
only.name    = 'Only defaults';
only.help    = {''};
only.strtype = 's';
only.num     = [1 Inf];
only.val     = {'yes'};

% exbranch only_defaults
only_defaults       = cfg_branch;
only_defaults.tag   = 'only_defaults';
only_defaults.name  = 'Load defaults only';
only_defaults.help  = {'Extracts and saves bval and bvec from the input NIfTI images. If the input is a set of 3D NIfTI images, it converts them into a single 4D NIfTI image (set filename A BIDS compliant filename can be set for the output 4D file.'};
only_defaults.val   = {only};


dummy_choice        = cfg_choice;
dummy_choice.tag    = 'dummy_choice';
dummy_choice.name   = 'Input type selection: Select existing data as input or extract bval/bvec and convert 3D NIfTIs';
dummy_choice.help   = {''};
dummy_choice.values = {extract, existing, existing_only_image, only_defaults};
% dummy_choice.val    = {extract};






% % ---------------------------------------------------------------------
% % Configure the acid Toolbox - load local, user-defined defaults file and
% % overwrite standard defaults
% % ---------------------------------------------------------------------
% acid_config         = cfg_exbranch;
% acid_config.tag     = 'acid_config';
% acid_config.name    = 'Configure toolbox';
% acid_config.val     = { acid_setdef };
% acid_config.help    = {['Customised default parameters can be set here by selecting ' ...
%     'a customised [acid_local_defaults_*.m] file. Type [help acid_local_defaults] for ' ...
%     'more details.']};
% acid_config.prog    = @acid_run_config;



%% satrtup
startup         = cfg_exbranch;
startup.tag     = 'startup';
startup.name    = 'Startup';
startup.val     = {dummy_choice, p_out, acid_setdef};
startup.help    = { '' };
startup.prog = @local_extract;
startup.vout = @vout_extract;


end

function out = local_extract(job)

tic

if isfield(job.dummy_choice,'existing') && isfield(job.dummy_choice.existing,'images_existing')
    job.filename = char(job.dummy_choice.existing.filename(1,:));
elseif isfield(job.dummy_choice,'extract') && isfield(job.dummy_choice.extract,'images')
    job.filename = char(job.dummy_choice.extract.filename(1,:));
elseif isfield(job.dummy_choice,'existing_only_image') && isfield(job.dummy_choice.existing_only_image,'images_existing')
    job.filename = char(job.dummy_choice.existing_only_image.filename(1,:));
elseif isfield(job.dummy_choice,'only_defaults')
else
    error('No filename was set!')
end

if isempty(job.p_out{1}) && ~isfield(job.dummy_choice,'only_defaults')
    if isfield(job.dummy_choice,'existing') && isfield(job.dummy_choice.existing,'images_existing')
        [p_out,~,~] = spm_fileparts(char(job.dummy_choice.existing.images_existing(1,:)));
    elseif isfield(job.dummy_choice,'extract') && isfield(job.dummy_choice.extract,'images')
        [p_out,~,~] = spm_fileparts(char(job.dummy_choice.extract.images(1,:)));
    elseif isfield(job.dummy_choice,'existing_only_image') && isfield(job.dummy_choice.existing_only_image,'images_existing')
        [p_out,~,~] = spm_fileparts(char(job.dummy_choice.existing_only_image.images_existing(1,:)));
    else
        error('No NIfTI file was selected!')
    end

    if isempty(job.filename)
        diary_path = [p_out filesep 'logfile_startup.txt'];
    else
        diary_path = [p_out filesep job.filename '-logfile_startup.txt'];
    end
    % start logging
    diary(diary_path)
    acid_startup_write_git_commitHash(p_out);
elseif isfield(job.dummy_choice,'only_defaults')


else
    p_out = job.p_out{1};
    if isempty(job.filename)
        diary_path = [p_out filesep 'logfile_startup.txt'];
    else
        diary_path = [p_out filesep job.filename '-logfile_startup.txt'];
    end

    % start logging
    diary(diary_path)
    acid_startup_write_git_commitHash(p_out);
end



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
    fprintf(['\nDefaults from "' char(deffnam{1}) '" are used.\n'])
else
    fprintf('\nStandard defaults are used.\n')
end
acid_get_defaults;


% read in bvals and bvecs
if isfield(job.dummy_choice,'existing')
    if isfield(job.dummy_choice.existing.bvals_bvecs,'bvals_bvecs_file')

        job.files1 = job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_file{1,1};


        try job.files2 = job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_file{2,1};
        catch
            job.files2 = job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_file{1,1};
        end



        [~,~,e] = fileparts(job.files1);

        if ~isempty(job.files1)
            switch lower(e)
                case '.mat'
                    if size(cell2mat(struct2cell(load(job.files1))),1) == 1 && size(cell2mat(struct2cell(load(job.files2))),1) == 3
                        job.bvals = cell2mat(struct2cell(load(job.files1)));
                        job.bvecs = cell2mat(struct2cell(load(job.files2)));
                    elseif size(cell2mat(struct2cell(load(job.files1))),1) == 3 && size(cell2mat(struct2cell(load(job.files2))),1) == 1
                        job.bvecs = cell2mat(struct2cell(load(job.files1)));
                        job.bvals = cell2mat(struct2cell(load(job.files2)));
                    elseif size(cell2mat(struct2cell(load(job.files1(1,:)))),1) == 4 || size(cell2mat(struct2cell(load(job.files1(1,:)))),2) == 4
                        bval_bvec_dummy = cell2mat(struct2cell(load(job.files1(1,:))));
                        job.bvecs = bval_bvec_dummy(2:4,:);
                        job.bvals = bval_bvec_dummy(1,:);
                    else
                        error('Unexpected file extension!')
                    end

                case '.nii'
                    error('Unexpected file extension!')
                otherwise
                    if size(dlmread(job.files1),1) == 1 && size(dlmread(job.files2),1) == 3
                        job.bvals = dlmread(job.files1);
                        job.bvecs = dlmread(job.files2);
                    elseif size(dlmread(job.files1),1) == 3 && size(dlmread(job.files2),1) == 1
                        job.bvecs = dlmread(job.files1);
                        job.bvals = dlmread(job.files2);
                    else
                        error('Unexpected file extension!')
                    end
            end

            P = char(job.dummy_choice.existing.images_existing);
            P = acid_load_4Dimage(P);
            filename = job.filename;
            bval = job.bvals;
            bvec = job.bvecs;
            normgrad = sqrt(sum(bvec .* bvec, 1));
            tolgrad    = 0.01; % tolerance for normalizing gradient directions
            MSKgrad = find(normgrad > 0);

            if(~isempty(find(normgrad(MSKgrad) ~= 1,1)))
                if(~isempty(find(abs(normgrad(MSKgrad)-1) > tolgrad,1))) % SM: if norm of gradients is within the tolerance message won't be visible / BF: And no changes on the b-vectors
                    warning('Diffusion gradients are now normalised!');
                    bvec(:, MSKgrad) = bsxfun(@rdivide, bvec(:, MSKgrad), normgrad(MSKgrad));
                end
            end
            clear normgrad;
            if size(bval,2) == size(bvec,2)
                spm_file_merge(P,[p_out filesep filename '.nii']);
                save([p_out filesep filename '_bval-bvec.mat'],'bval','bvec')
                Mtmp1 = [p_out filesep filename '.bval'];
                Mtmp2 = [p_out filesep filename '.bvec'];
                dlmwrite(Mtmp1, bval);
                dlmwrite(Mtmp2, bvec);
            else
                error('Number of b-vectors and b-values are not the same!')
            end


        else
            job.bvals = '';
            job.bvecs = '';
        end

    elseif isfield(job.dummy_choice.existing.bvals_bvecs,'bvals_bvecs_exp_type')

        if isfield(job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_exp_type,'bvals_exp')
            if ~isempty(job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_exp_type.bvals_exp)
                job.bvals = job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_exp_type.bvals_exp;
            else
                job.bvals = '';
            end
        end

        if isfield(job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_exp_type,'bvecs_exp')
            if ~isempty(job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_exp_type.bvecs_exp)
                job.bvecs = job.dummy_choice.existing.bvals_bvecs.bvals_bvecs_exp_type.bvecs_exp;
            else
                job.bvecs = '';
            end
        end


        P = char(job.dummy_choice.existing.images_existing);
        P = acid_load_4Dimage(P);
        filename = job.filename;
        bval = job.bvals;
        bvec = job.bvecs;

        normgrad = sqrt(sum(bvec .* bvec, 1));
        tolgrad    = 0.01; % tolerance for normalizing gradient directions
        MSKgrad = find(normgrad > 0);

        if(~isempty(find(normgrad(MSKgrad) ~= 1,1)))
            if(~isempty(find(abs(normgrad(MSKgrad)-1) > tolgrad,1))) % SM: if norm of gradients is within the tolerance message won't be visible / BF: And no changes on the b-vectors
                warning('Diffusion gradients are now normalised!');
                bvec(:, MSKgrad) = bsxfun(@rdivide, bvec(:, MSKgrad), normgrad(MSKgrad));
            end
        end
        clear normgrad;

        if size(bval,2) == size(bvec,2)
            spm_file_merge(P,[p_out filesep filename '.nii']);
            save([p_out filesep filename '_bval-bvec.mat'],'bval','bvec')
            Mtmp1 = [p_out filesep filename '.bval'];
            Mtmp2 = [p_out filesep filename '.bvec'];
            dlmwrite(Mtmp1, bval);
            dlmwrite(Mtmp2, bvec);
        else
            error('Number of b-vectors and b-values are not the same!')
        end
    end
elseif isfield(job.dummy_choice,'extract')

    [bval, bvec] = acid_extract_bval_bvec(char(job.dummy_choice.extract.images), job.filename, p_out);

    job.bvals = bval;
    job.bvecs = bvec;

elseif isfield(job.dummy_choice,'only_defaults')

    job.bvals    = '';
    job.bvecs    = '';
    job.filename = '123';
    p_out        = '';

elseif isfield(job.dummy_choice,'existing_only_image')

    P = char(job.dummy_choice.existing_only_image.images_existing);
    P = acid_load_4Dimage(P);
    filename = job.filename;
    job.bvals    = '';
    job.bvecs    = '';
    spm_file_merge(P,[p_out filesep filename '.nii']);

end


out.rsource{1} = [p_out filesep job.filename '.nii,1'];

out.bvals = job.bvals;
out.bvecs = job.bvecs;

T = toc/60;

T = duration(minutes(T),'format','hh:mm:ss');
disp(['The total time for startup was: ' char(T) '.']);
diary off

end

function dep = vout_extract(job)

dummy_bval_bvec = 1;
dummy_all = 1;
dummy_extracted = 0;



if isfield(job.dummy_choice, 'existing_only_image')
    dummy_bval_bvec = 0;
elseif isfield(job.dummy_choice, 'only_defaults')
    dummy_all = 0;
elseif isfield(job.dummy_choice, 'extract')
    dummy_extracted = 1;
end

if dummy_all == 1
    dep(1)            = cfg_dep;
    if dummy_extracted == 1
        dep(1).sname      = 'Created 4D NIfTI file';
    else
        dep(1).sname      = 'Imported 4D NIfTI file';
    end
    dep(1).src_output = substruct('.','rsource');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if dummy_bval_bvec == 1
        dep(2)            = cfg_dep;
        if dummy_extracted == 1
            dep(2).sname      = 'Extracted b-values';
        else
            dep(2).sname      = 'Imported b-values';
        end
        dep(2).src_output = substruct('.','bvals');
        dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
        dep(3)            = cfg_dep;
        if dummy_extracted == 1
            dep(3).sname      = 'Extracted b-vectors';
        else
            dep(3).sname      = 'Imported b-vectors';
        end
        dep(3).src_output = substruct('.','bvecs');
        dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});
    end
else
    dep(1)            = cfg_dep;
    dep(1).sname      = 'No dependency!';
    dep(1).src_output = substruct('.','');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
end
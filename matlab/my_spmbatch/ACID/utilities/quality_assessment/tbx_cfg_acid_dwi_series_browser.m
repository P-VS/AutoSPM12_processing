function browser = tbx_cfg_acid_dwi_series_browser(bvals_bvecs)

% source image
sources         = cfg_files;
sources.tag     = 'sources';
sources.name    = 'dMRI dataset';
sources.help    = {'Select the dMRI dataset.'};
sources.filter  = 'image';
sources.ufilter = '.*';
sources.num     = [1 Inf];

% exbranch:
browser      = cfg_exbranch;
browser.tag  = 'browser';
browser.name = 'DWI series browser';
browser.help = {'dMRI dataset browser.'};
browser.val  = {sources, bvals_bvecs};
browser.prog = @local_browser;
browser.vout = @out_browser;

end

function out = local_browser(job)

    % read in bvals and bvecs
    if isfield(job.bvals_bvecs,'bvals_bvecs_file')

        job.sources = cell2mat(job.bvals_bvecs.bvals_bvecs_file);
        [~,~,e] = fileparts(job.sources(1,:));

        if ~isempty(job.sources(1,:))
            switch lower(e)
                case '.mat'
                    if size(cell2mat(struct2cell(load(job.sources(1,:)))),1) == 1 || size(cell2mat(struct2cell(load(job.sources(1,:)))),2) == 1
                        job.bvals = cell2mat(struct2cell(load(job.sources(1,:))));
                        job.bvecs = cell2mat(struct2cell(load(job.sources(2,:))));
                    elseif size(cell2mat(struct2cell(load(job.sources(1,:)))),1) == 3 || size(cell2mat(struct2cell(load(job.sources(1,:)))),2) == 3
                        job.bvecs = cell2mat(struct2cell(load(job.sources(1,:))));
                        job.bvals = cell2mat(struct2cell(load(job.sources(2,:))));
                    else
                        error('Unexpected file extension!')
                    end

                case '.nii'
                    error('Unexpected file extension!')
                otherwise
                    if size(dlmread(job.sources(1,:)),1) == 1 || size(dlmread(job.sources(1,:)),2) == 1
                        job.bvals = dlmread(job.sources(1,:));
                        job.bvecs = dlmread(job.sources(2,:));
                    elseif size(dlmread(job.sources(1,:)),1) == 3 || size(dlmread(job.sources(1,:)),2) == 3
                        job.bvecs = dlmread(job.sources(1,:));
                        job.bvals = dlmread(job.sources(2,:));
                    else
                        error('Unexpected file extension!')
                    end
            end
        else
            job.bvals = '';
            job.bvecs = '';
        end

    elseif isfield(job.bvals_bvecs,'bvals_bvecs_exp_type')

        if isfield(job.bvals_bvecs.bvals_bvecs_exp_type,'bvals_exp')
            if ~isempty(job.bvals_bvecs.bvals_bvecs_exp_type.bvals_exp)
                job.bvals = job.bvals_bvecs.bvals_bvecs_exp_type.bvals_exp;
            else
                job.bvals = '';
            end
        end

        if isfield(job.bvals_bvecs.bvals_bvecs_exp_type,'bvecs_exp')
            if ~isempty(job.bvals_bvecs.bvals_bvecs_exp_type.bvecs_exp)
                job.bvecs = job.bvals_bvecs.bvals_bvecs_exp_type.bvecs_exp;
            else
                job.bvecs = '';
            end
        end
    end

    acid_dwi_series_browser(char(job.sources), job.bvals, job.bvecs);

    out.files   = job.sources;
    out.matfile = acid_spm_file(char(job.sources(1)),'prefix','labeled_slices_','format','.mat');
    if isfield(job.bvals_bvecs,'bvals_bvecs_file')
        out.bval = cellstr(job.bvals_type.browser_bvals_file);
        out.bvec = cellstr(job.bvecs_type.bvecs_file);
    elseif isfield(job.bvals_bvecs,'bvals_bvecs_exp')
        out.bval = job.bvals_type.bvals_exp;
        out.bvec = job.bvecs_type.bvecs_exp;
    end
end

function dep = out_browser(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Images';
    dep(1).src_output = substruct('.','files');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(2)            = cfg_dep;
    dep(2).sname      = 'b-values';
    dep(2).src_output = substruct('.','bval');
    dep(2).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    dep(3)            = cfg_dep;
    dep(3).sname      = 'b-vectors';
    dep(3).src_output = substruct('.','bvec');
    dep(3).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    dep(4)            = cfg_dep;
    dep(4).sname      = 'Labeled volumes';
    dep(4).src_output = substruct('.','matfile');
    dep(4).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
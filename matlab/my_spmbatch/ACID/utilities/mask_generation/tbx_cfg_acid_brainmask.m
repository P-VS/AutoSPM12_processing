function make_brainMSK = tbx_cfg_acid_brainmask

% image to segment
P         = cfg_files;
P.tag     = 'P';
P.name    = 'Brain image (.nii)';
P.help    = {'Select a image for segmentation and mask generation. MD images can not be used because the voxel values are to low fors SPMs segmentation function.'};
P.filter  = 'image';
P.ufilter = '.*';
P.num     = [0 1];

% tissue probability maps (TPMs) of GM, WM, and CSF
TPMs         = cfg_files;
TPMs.tag     = 'TMPs';
TPMs.name    = 'c1-c3 segmentations';
TPMs.help    = {'...'};
TPMs.filter  = 'image';
TPMs.ufilter = '.*';
TPMs.num     = [0 3];

% input type for brain mask generation
input_choice        = cfg_choice;
input_choice.tag    = 'input_choice';
input_choice.name   = 'Select input type for brain mask generation';
input_choice.help   = {''};
input_choice.values = {P, TPMs};
input_choice.val    = {P};

% images to be masked
P_others         = cfg_files;
P_others.tag     = 'P_others';
P_others.name    = 'DTI maps to be masked';
P_others.help    = {'...'};
P_others.filter  = 'image';
P_others.ufilter = '.*';
P_others.num     = [0 Inf];
P_others.val     = {{''}};

% Inerpolation order
perc         = cfg_entry;
perc.tag     = 'perc';
perc.name    = 'Brain coverage';
perc.help    = {'Brain coverage (0: brain covers full field-of-view).'};
perc.strtype = 'e';
perc.num     = [1 1];
perc.val     = {0.8};

% smoothing kernel used in both, DTI fit and brain mask
smk         = cfg_entry;
smk.tag     = 'smk';
smk.name    = 'Smoothing kernel';
smk.help    = {'This option reduces holes in the brain mask by smoothing the mask or the tissue probability maps. A value of zero indicates no smoothing.'};
smk.strtype = 'e';
smk.num     = [1 3];
smk.val     = {[3 3 3]};

% exbranch
make_brainMSK      = cfg_exbranch;
make_brainMSK.tag  = 'make_brainMSK';
make_brainMSK.name = 'Create brain mask';
make_brainMSK.val  = {input_choice, P_others, perc, smk};
make_brainMSK.help = {
    'This function constructs a mask from an input image volume or c1-c3 segmentations and applies the mask on other images. Recommendation about image(s) that are used for mask creation: apply ???SPM New Segment??? on low-b-value image and use take c1-c3 for mask creation.'
    }';
make_brainMSK.prog = @local_make_brainMSK;
make_brainMSK.vout = @vout_make_brainMSK;

end

function out = local_make_brainMSK(job)

    dummy_options.smk  = job.smk;
    dummy_options.perc = job.perc;

    if isfield(job.input_choice,'P')
        dummy_segmentation = 1;
        job.P = job.input_choice.P;
    else
        dummy_segmentation = 0;
        job.P = job.input_choice.TPMs;
    end

    Ni = acid_brainmask(char(job.P), char(job.P_others), [], dummy_options, dummy_segmentation);
    out.mask = {[Ni.dat.fname ',' num2str(1)]};

    try dummy_depend = ~isempty(job.P_others{1,1});
        if dummy_depend == 1
            
            keyword = 'BRAIN-MASKED';
            p = fileparts(Ni.dat.fname);
            
            for i = 1:size(job.P_others,1)
                fname = acid_bids_filename(spm_vol(char(job.P_others(i,:))), keyword, 'same', '.nii');
                tmp{i,1} = [p filesep fname ',1'];
            end
            out.masked = tmp(:);
        end
    catch
    end
end

function dep = vout_make_brainMSK(job)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Brain mask';
    dep(1).src_output = substruct('.','mask');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        
    try dummy_depend = ~isempty(job.P_others{1,1});     
        if dummy_depend == 1
            dep(2)            = cfg_dep;
            dep(2).sname      = 'Brain masked images';
            dep(2).src_output = substruct('.','masked');
            dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        end
    catch
    end 
end
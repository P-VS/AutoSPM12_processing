function [dummy_Tfreiburg, mask, bvals_bvecs, in_vols, RMatrix, p_out, bvals] = tbx_cfg_acid_params

%% in_vols
in_vols         = cfg_files;
in_vols.tag     = 'in_vols';
in_vols.name    = 'Input images';
in_vols.help    = {'Select the dMRI data.'};
in_vols.filter  = 'image';
in_vols.ufilter = '.*';
in_vols.num     = [0 Inf];

%% b values / b vectors
bvals_bvecs_file         = cfg_files;
bvals_bvecs_file.tag     = 'bvals_bvecs_file';
bvals_bvecs_file.name    = 'Select b-values and b-vectors files (*.mat|txt|bval|bvec)';
bvals_bvecs_file.help    = {'Select two files specifying the b-values (.bval) and b-vectors (.bvec) of the N dMRI volumes.'
    'Note: The file should contain a 1 x N vector, where b-values should appear in the same order as the source images were entered. b-values is expected in units of s/mm^2.'};
bvals_bvecs_file.filter  = 'any';
bvals_bvecs_file.ufilter = 'mat|txt|bval|bvec';
bvals_bvecs_file.num     = [0 2];

bvals_file         = cfg_files;
bvals_file.tag     = 'bvals_file';
bvals_file.name    = 'Select b-values file (*.mat|txt|bval)';
bvals_file.help    = {'Select two files specifying the b-values (.bval) and b-vectors (.bvec) of the N dMRI volumes.'
    'Note: The file should contain a 1 x N vector, where b-values should appear in the same order as the source images were entered. b-values is expected in units of s/mm^2.'};
bvals_file.filter  = 'any';
bvals_file.ufilter = 'mat|txt|bval';
bvals_file.num     = [0 1];

bvals_exp         = cfg_entry;
bvals_exp.tag     = 'bvals_exp';
bvals_exp.name    = 'b-values';
bvals_exp.help    = {'Specify the b-values (bval) of the N dMRI volumes.'
    'Note: The expression should contain a 1 x N vector, where b-values should appear in the same order as the source images were entered (b-values: 1 x N, b-vectors: (2-4) x N). b-values are expected in units of s/mm^2.'};
bvals_exp.strtype = 'e';
bvals_exp.num     = [Inf Inf];

bvecs_exp         = cfg_entry;
bvecs_exp.tag     = 'bvecs_exp';
bvecs_exp.name    = 'b-vectors';
bvecs_exp.help    = {'Specify the b-vectors (bvec) of the N dMRI volumes.'
    'Note: The expression should contain a 3 x N vector, where b-vectors should appear in the same order as the source images were entered. b-vectors are expected in units of s/mm^2.'};
bvecs_exp.strtype = 'e';
bvecs_exp.num     = [Inf Inf];

bvals_bvecs_exp_type         = cfg_branch;
bvals_bvecs_exp_type.tag     = 'bvals_bvecs_exp_type';
bvals_bvecs_exp_type.name    = 'Expression/Dependency';
% bvals_bvecs_exp_type.help    = {'Specify the b-values (bval) and b-vectors of the N dMRI volumes.'
%     'Note: The expression should contain a 4 x N vector, where b-values should appear in the same order as the source images were entered (b-values: 1 x N, b-vectors: (2-4) x N). b-values are expected in units of s/mm^2.'};
bvals_bvecs_exp_type.val  = {bvals_exp, bvecs_exp};

bvals_bvecs        = cfg_choice;
bvals_bvecs.tag    = 'bvals_bvecs';
bvals_bvecs.name   = 'Select b-values and b-vectors input type';
bvals_bvecs.help   = {''};
bvals_bvecs.values = {bvals_bvecs_file, bvals_bvecs_exp_type};
bvals_bvecs.val    = {bvals_bvecs_exp_type};

bvals        = cfg_choice;
bvals.tag    = 'bvals';
bvals.name   = 'Select b-values input type';
bvals.help   = {''};
bvals.values = {bvals_file, bvals_exp};
bvals.val    = {bvals_exp};

%% b values
% b_vals         = cfg_entry;
% b_vals.tag     = 'b_vals';
% b_vals.name    = 'b-values (bval)';
% b_vals.help    = {'Provide an 1 x N  - array with b-values, b-values should appear in the same order as the low- and high-diffusion weighted images were entered. b-values is expected in units of s/mm^2.' 
%                   'Note that the provided entry is only for illustration.'};
% b_vals.strtype = 'e';
% b_vals.num     = [Inf Inf];
% b_vals.val     = {[5 1000 1000 2000]};

%% Reorientation Matrix for b-vectors
RMatrix         = cfg_entry;
RMatrix.tag     = 'RMatrix';
RMatrix.name    = 'Reorientation Matrix';
RMatrix.help    = {
                      'If the vendor uses another coordinate system than the coordinate system, in which your b-vectors were defined, you need to reorient them.'
                      'Provide a 3 x 3  - matrix to reorient b-vectors.'
};
RMatrix.strtype = 'e';
RMatrix.num     = [3 3];
RMatrix.val     = {[1 0 0; 0 1 0; 0 0 1]};

%% Binary Map defining Region of Interest
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Region of interest';
mask.val     = {{''}};
mask.help    = {'The computeation is restricted to the binary region of interest mask. It is recommened providing a mask for large datasets to avoid memory issues.'};
mask.filter  = 'image';
mask.ufilter = '.*';
mask.num     = [0 Inf];

%% dummy_robust
dummy_Tfreiburg   = cfg_menu;
dummy_Tfreiburg.tag     = 'dummy_Tfreiburg';
dummy_Tfreiburg.name    = 'Write Freiburg Tractography format';
dummy_Tfreiburg.help    = {'Write Freiburg Tractography format. The path of the Freiburg Tractography tools have to be included - otherwise it will not work.'};
dummy_Tfreiburg.labels = {
               'NO'
               'YES'
}';
dummy_Tfreiburg.values = {0 1};
dummy_Tfreiburg.val    = {0};

%% output directory
p_out         = cfg_files;
p_out.tag     = 'p_out';
p_out.name    = 'Output directory';
p_out.help    = {'Select the output directory. If unspecified, the output will be saved in the derivatives folder.'};
p_out.filter  = 'dir';
p_out.ufilter = '.*';
p_out.num     = [0 1];
p_out.val     = {{''}};

end
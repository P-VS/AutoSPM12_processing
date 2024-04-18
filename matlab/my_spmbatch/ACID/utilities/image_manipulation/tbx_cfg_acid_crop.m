function crop_new = tbx_cfg_acid_crop

% images to crop
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Input images';
images.help    = {'Select the dMRI dataset you want to crop. You can load in either a set of 3d nifti or a single 4d nifti file.'};
images.filter  = 'image';
images.ufilter = '.*';
images.num     = [0 Inf];

% new matrix size
crop_new_matrix         = cfg_entry;
crop_new_matrix.tag     = 'crop_new_matrix';
crop_new_matrix.name    = 'New image matrix';
crop_new_matrix.help    = {'Enter the matrix size of the cropped image.'
    'The three values represent matrix size in the x, y, and z dimensions, respectively'};
crop_new_matrix.strtype = 'e';
crop_new_matrix.num     = [1 3];
crop_new_matrix.val     = {};

% center of cropping: x and y voxel coordinates
crop_new_midxy         = cfg_entry;
crop_new_midxy.tag     = 'crop_new_midxy';
crop_new_midxy.name    = 'Center of cropping: x and y voxel coordinates';
crop_new_midxy.help    = {'Enter the x and y voxel coordinates around which the image(s) will be cropped.'
    'If left unspecified, the midpoint has to be selected manually in a pop-up window.'};
crop_new_midxy.strtype = 'e';
crop_new_midxy.num     = [Inf Inf];
crop_new_midxy.val     = {[]};

% center of cropping: z voxel coordinate
crop_new_midz         = cfg_entry;
crop_new_midz.tag     = 'crop_new_midz';
crop_new_midz.name    = 'Center of cropping: z voxel coordinate';
crop_new_midz.help    = {'Enter the z voxel coordinate around which the image(s) will be cropped.'};
crop_new_midz.strtype = 'e';
crop_new_midz.num     = [Inf Inf];
crop_new_midz.val     = {[]};

% exbranch
crop_new      = cfg_exbranch;
crop_new.tag  = 'crop_new';
crop_new.name = 'Cropping: create new';
crop_new.help = {'Create new cropping'};
crop_new.val  = {images, crop_new_matrix, crop_new_midxy, crop_new_midz};
crop_new.prog = @local_crop_images;
crop_new.vout = @vout_crop_images;

end

function out = local_crop_images(job)
    [Vout, fname_params] = acid_crop(char(job.images), job.crop_new_matrix, job.crop_new_midxy, job.crop_new_midz);

    out.cfiles  = {[Vout(1).fname ',1']};  
    out.matfile = {fname_params};

end

function dep = vout_crop_images(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Cropped images';
    dep(1).src_output = substruct('.','cfiles');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(2)            = cfg_dep;
    dep(2).sname      = 'Cropping parameters';
    dep(2).src_output = substruct('.','matfile');
    dep(2).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end
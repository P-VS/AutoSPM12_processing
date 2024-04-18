function crop_apply = tbx_cfg_acid_crop_apply

% images to crop
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Images to crop';
images.help    = {'Select the dMRI dataset you want to crop. You can load in either a set of 3D nifti or a single 4D nifti file.'};
images.filter  = 'image';
images.ufilter = '.*';
images.num     = [0 Inf];

% cropping parameters
params         = cfg_files;
params.tag     = 'params';
params.name    = 'Cropping parameters';
params.help    = {'Select the matfile containing the cropping parameters'};
params.filter  = 'mat';
params.ufilter = '.*';
params.num     = [0 Inf];

% exbranch
crop_apply      = cfg_exbranch;
crop_apply.tag  = 'crop_apply';
crop_apply.name = 'Cropping: apply existing';
crop_apply.help = {'Apply existing cropping'};
crop_apply.val  = {images, params};
crop_apply.prog = @local_crop_images;
crop_apply.vout = @vout_crop_images;
end

function out = local_crop_images(job)

    params = load(char(job.params));
    Vout = acid_crop(char(job.images), params.matrix, params.mid_xy, params.mid_z);
    
    for i = 1:size(Vout,1)
        [a,b,c] = fileparts(Vout(1).fname);  
        fname = [b];    
        out.cfiles{i,:} = [a filesep fname c ',' num2str(i)]; 
    end

    out.matrix = params.matrix;
    out.mid_z  = params.mid_z;
    out.mid_xy = params.mid_xy;
end

function dep = vout_crop_images(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Cropped images';
    dep(1).src_output = substruct('.','cfiles');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(2)            = cfg_dep;
    dep(2).sname      = 'New image matrix';
    dep(2).src_output = substruct('.','matrix');
    dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
    dep(3)            = cfg_dep;
    dep(3).sname      = 'Middle slice';
    dep(3).src_output = substruct('.','mid_z');
    dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});  
    dep(4)            = cfg_dep;
    dep(4).sname      = 'Midpoint in the plain';
    dep(4).src_output = substruct('.','mid_xy');
    dep(4).tgt_spec   = cfg_findspec({{'strtype','e'}});
end
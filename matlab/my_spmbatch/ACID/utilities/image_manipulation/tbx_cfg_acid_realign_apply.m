function realign_apply = tbx_cfg_acid_realign_apply

% target
target         = cfg_files;
target.tag     = 'target';
target.name    = 'Target images';
target.help    = {'The target image is the image to which the source image is registered. The number of slicewise transformations will be based on the dimension of this image.'};
target.filter  = 'image';
target.ufilter = '.*';
target.num     = [1 1];

% other images
images         = cfg_files;
images.tag     = 'images';
images.name    = 'Source images';
images.help    = {'Select images to apply slicwise transformation (should be alligned with the source image).'};
images.filter  = 'image';
images.ufilter = '.*';
images.num     = [1 Inf];

% parameters
params         = cfg_files;
params.tag     = 'params';
params.name    = 'Matfile of slicewise transformations';
params.help    = {'Select matfile(s) that was/were produced when pressing "apply" on manuel registration. Note that you can have more than one matfile.'};
params.filter  = 'mat';
params.ufilter = '.*';
params.num     = [1 Inf];

% exbranch
realign_apply         = cfg_exbranch;
realign_apply.tag     = 'realign_apply';
realign_apply.name    = 'Slice-wise realign: apply existing';
realign_apply.val     = {target, images, params};
realign_apply.help    = {
                    ''
};
realign_apply.prog    = @local_apply_realign_slices;
realign_apply.vout    = @vout_apply_realign_slices;

end

function out = local_apply_realign_slices(job)

    Vout = acid_realign_apply(char(job.target), char(job.images), char(job.params)); 

    for i = 1:size(Vout,1)
        [a,b,c] = fileparts(Vout(1).fname);
        fname = b;
        out.rfiles{i,:} = [a filesep fname c ',' num2str(i)];
    end

end

function dep = vout_apply_realign_slices(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Realigned images';
    dep(1).src_output = substruct('.','rfiles');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    
end
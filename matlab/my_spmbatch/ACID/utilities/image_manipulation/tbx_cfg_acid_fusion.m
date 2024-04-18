function fusion = tbx_cfg_acid_fusion

% source image
target         = cfg_files;
target.tag     = 'target';
target.name    = 'Target image';
target.help    = {''};
target.filter  = 'image';
target.ufilter = '.*';
target.num     = [0 1];

% brain image
b_image         = cfg_files;
b_image.tag     = 'b_image';
b_image.name    = 'Brain FA map';
b_image.help    = {''};
b_image.filter  = 'image';
b_image.ufilter = '.*';
b_image.num     = [0 1];

% spinal cord image
sc_image         = cfg_files;
sc_image.tag     = 'sc_image';
sc_image.name    = 'Spinal-Cord FA map';
sc_image.help    = {''};
sc_image.filter  = 'image';
sc_image.ufilter = '.*';
sc_image.num     = [0 1];

% % eigenvector images
% eigenvector_brain_image         = cfg_files;
% eigenvector_brain_image.tag     = 'eigenvector_brain_image';
% eigenvector_brain_image.name    = 'Brain eigenvector images (V1,V2,V3 order is mandatory!)';
% eigenvector_brain_image.help    = {''};
% eigenvector_brain_image.filter  = 'image';
% eigenvector_brain_image.ufilter = '.*';
% eigenvector_brain_image.num     = [3 3];
% 
% % eigenvector images
% eigenvector_sc_image         = cfg_files;
% eigenvector_sc_image.tag     = 'eigenvector_sc_image';
% eigenvector_sc_image.name    = 'Spinal cord eigenvector images (V1,V2,V3 order is mandatory!)';
% eigenvector_sc_image.help    = {''};
% eigenvector_sc_image.filter  = 'image';
% eigenvector_sc_image.ufilter = '.*';
% eigenvector_sc_image.num     = [3 3];
    
% voxel dimensions
overlap         = cfg_entry;
overlap.tag     = 'overlap';
overlap.name    = 'Overlap correction';
overlap.help    = {''};
overlap.strtype = 'e';
overlap.num     = [1 1];
overlap.val     = {3};
    
% interpolation
interp         = cfg_entry;
interp.tag     = 'interp';
interp.name    = 'Interpolation order';
interp.help    = {'Interpolation order as defined in spm_reslice. A splin-interpolation is used by default.'};
interp.strtype = 'e';
interp.num     = [1 1];
interp.val     = {-7};

%% Exbranch
fusion       = cfg_exbranch;
fusion.tag   = 'fusion';
fusion.name  = 'Fusion';
fusion.val   = {target b_image sc_image overlap interp};
fusion.help  = {'...'};
                    
fusion.prog  = @local_fusion;
% fusion.vout  = @vout_fusion;

end

function local_fusion(job)



acid_fusion(job.target, job.b_image, job.sc_image, job.overlap, job.interp);


end
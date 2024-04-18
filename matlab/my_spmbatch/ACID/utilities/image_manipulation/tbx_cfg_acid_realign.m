function realign = tbx_cfg_acid_realign

% source image
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source image';
source.help    = {'Select the image to be realigned.'};
source.filter  = 'image';
source.ufilter = '.*';
source.num     = [1 1];

% target image
target         = cfg_files;
target.tag     = 'target';
target.name    = 'Target image';
target.help    = {'Select the reference image. This image will remain stationary during the realignment process.'};
target.filter  = 'image';
target.ufilter = '.*';
target.num     = [1 1];

% other images
others = cfg_files;
others.tag     = 'others';
others.name    = 'Other images';
others.help    = {'Select images to apply slice-wise transformation (should be alligned with the source image).'};
others.filter  = 'image';
others.ufilter = '.*';
others.num     = [0 Inf];

% exbranch
realign = cfg_exbranch;
realign.tag    = 'realign';
realign.name   = 'Slice-wise realign: create new';
realign.help   = {'Create new slice-wise realignment'};
realign.val = {source, target, others};
realign.prog = @local_realign_slices;
realign.vout = @vout_realign_slices;

end

function out = local_realign_slices(job)
    
    p_out = acid_realign(char(job.source), char(job.target), char(job.others));

    keyword = 'REALIGN';
    
    % get output parameter file
    fname_params = acid_bids_filename(spm_vol(char(job.source)),[keyword '-params'],'','.mat');
    out.params = [p_out filesep fname_params];
    
    % get output realigned image(s)
    P = char(job.source);
    fname = acid_bids_filename(spm_vol(P),keyword,'same','.nii');
    out.rfiles{1,:} = [p_out filesep fname P(1,strfind(P,','):end)];
    
    if ~isempty(job.others)
        P = char(job.others);
        for i = 2:size(job.others,1)+1      
            fname = acid_bids_filename(spm_vol(P(i-1,:)),keyword,'same','.nii');
            out.rfiles{i,:} = [p_out filesep fname P(i-1,strfind(P(i-1,:),','):end)];
        end
        if strcmp(out.rfiles{1,:},out.rfiles{2,:})
            out.rfiles(1) = [];
        end
    end
    
end

function dep = vout_realign_slices(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Transformation parameters';
    dep(1).src_output = substruct('.','params');
    dep(1).tgt_spec   = cfg_findspec({{'filter','r','strtype','e'}});
    dep(2)            = cfg_dep;
    dep(2).sname      = 'Realigned images';
    dep(2).src_output = substruct('.','rfiles');
    dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
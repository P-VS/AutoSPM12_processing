function hysco_choice = tbx_cfg_acid_hysco(in_vols, dummy_hysco_apply_multi)

% reference blip-up
source_bu      = in_vols;
source_bu.tag  = 'source_bu';
source_bu.name = 'Reference blip-up image';
source_bu.help = {
''
'Select one image volume acquired with blip-up. The field inhomogeneity is estimated by minimizing the sum-of-squared difference between this image and the blip-down image chosen below and regularization.'
''};
 source_bu.num = [0 1];

% reference blip-down
source_bd      = in_vols;
source_bd.tag  = 'source_bd';
source_bd.name = 'Reference blip-down image';
source_bd.help = {''
'Select one image volume acquired with blip-down. The field inhomogeneity is estimated by minimizing the sum-of-squared difference between this image and the blip-up image chosen above and regularization.'
''};
source_bd.num  = [0 1];

% other blip-up image(s)
others_bu      = in_vols;
others_bu.tag  = 'others_bu';
others_bu.name = 'Other blip-up image(s)';
others_bu.help = {''
'(optional) Choose other image volumes acquired with blip-up that need to be corrected. The data is corrected by applying the transformation estimated by the reference bip-up/down data. If an equal number of blip-up and blip-down data is provided, you may also want to disable ''Apply to other images''.'
''};
others_bu.val  = {{''}};

% other blip-down image(s)
others_bd      = in_vols;
others_bd.tag  = 'others_bd';
others_bd.name = 'Other blip-down image(s)';
others_bd.help = {''
'(optional) Choose other image volumes acquired with blip-down that need to be corrected. The data is corrected by applying the transformation estimated by the reference bip-up/down data. If an equal number of blip-up and blip-down data is provided, you may also  want to disable ''Apply to other images''.'
''};
others_bd.val  = {{''}};

% phase-encoding direction
phdir        = cfg_menu;
phdir.tag    = 'phdir';
phdir.name   = 'Phase-encoding direction';
phdir.help   = {''
'Specify the phase-encoding direction.'
''};
phdir.labels = {'x','y','z'};
phdir.values = {1 2 3};
phdir.val    = {2};

% dummy_fast: maximal data resolution
dummy_fast        = cfg_menu;
dummy_fast.tag    = 'dummy_fast';
dummy_fast.name   = 'Maximal data resolution';
dummy_fast.help   = {''
'Choose the finest discretization level for field inhomogeneity estimation. If set to ''full'' a multi-level strategy with three discretization levels is used, where the resolution on the finest level equals the data resolution. To save computation time, choose ''half''. The multi-level scheme will be stopped after the second level (i.e. half of data resolution) and the inhomogeneity estimate will be interpolated to the data resolution.'
''};
dummy_fast.labels = {
               'half'
               'full'
};
dummy_fast.values = {0 1};
dummy_fast.val    = {1};

%% HySCO: apply existing
% field inhomogeneity map
fmap         = cfg_files;
fmap.tag     = 'fmap';
fmap.name    = 'Field inhomogeneity map';
fmap.help    = {
''
'Select the map of field inhomogeneity.'
''};
fmap.filter  = 'image';
fmap.ufilter = '.*';
fmap.num     = [0 Inf];

%% Exbranch
hysco      = cfg_exbranch;
hysco.tag  = 'hysco';
hysco.name = 'HySCO: create new';
hysco.val  = {source_bu source_bd others_bu others_bd phdir dummy_fast};
hysco.help = {
                    'Hyperelastic susceptibility artifact correction of diffusion weighted images'
                    ''
                    'HySCO enables correction of susceptibility artifacts in diffusion weighted images using two images acquired with reversed phase encoding gradients.'
                    ''
                    'A nonlinear regularization functional, which is inspired by hyperelasticity, ensures smoothness of the field inhomogeneity and invertibility of the geometrical transformations. It consists of two components for which regularization parameters have to be chosen: A "diffusion" part that enforces smoothness of the field and a "Jacobian" part that guarantees invertibility of the estimated transformations.'
                    ''
                    'This code is faster and more accurate than HySCO 1.0 due to the inexact Gauss-Newton optimization proposed in Macdonald and Ruthotto 2016'
                    ''
                    'For more information and references see:'
                    ''
                    'MacDonald, J, Ruthotto, L, '
                    'Efficient Numerical Optimization For Susceptibility Artifact Correction Of EPI-MRI'
                    'arXiv, 2016'
                    ''
                    'Ruthotto, L, Kugel, H, Olesch, J, Fischer, B, Modersitzki, J, Burger, M, and Wolters, C H. '
                    'Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images. '
                    'Physics in Medicine and Biology, 57(18), 5715-5731; 2012.'
                    ''
                    'Ruthotto, L, Mohammadi, S, Heck, C, Modersitzki, J, and Weiskopf, N.'
                    'HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM.'
                    'Presented at the Bildverarbeitung fuer die Medizin 2013.'
};
hysco.prog = @local_hysco;
hysco.vout = @vout_hysco;

%% HySCO: apply existing 

% define space defining image 
source_bu_apply      = source_bu;
source_bu_apply.help = {
    'Space defining image. Select the image volume acquired with blip-up that was used as reference image to estimate fieldmap.'};

% hysco_apply: apply susceptibility correction
hysco_apply      = cfg_exbranch;
hysco_apply.tag  = 'hysco_apply';
hysco_apply.name = 'HySCO: apply existing';
hysco_apply.val  = {source_bu_apply others_bu others_bd fmap phdir};
hysco_apply.help = {
                    'This option applies the field-inhomogeneities from the Hyperelastic susceptibility artifact correction of diffusion weighted images to undistort images of choice.'					
};
hysco_apply.prog = @local_hysco_apply;
hysco_apply.vout = @vout_hysco_apply;


%% HySCO: ???
hysco_apply_multi      = cfg_exbranch;
hysco_apply_multi.tag  = 'hysco_apply_multi';
hysco_apply_multi.name = 'HySCO: Write corrected images using multiple field maps';
hysco_apply_multi.val  = {source_bu_apply others_bu others_bd fmap phdir};
hysco_apply_multi.help = {
                    'This module writes the Jacobian estimated from the HySCO map, as well as unwarped blip-up and down image.'
                    'If more than one HySCO map is provided, the HySCO maps will be added up.'
};
hysco_apply_multi.prog = @local_hysco_write_multi;
hysco_apply_multi.vout = @vout_hysco_write_multi;

%% HySCO: 

% unwarped blip-up images
hysco_bu         = cfg_files;
hysco_bu.tag     = 'hysco_bu';
hysco_bu.name    = 'Unwarped blip-up images';
hysco_bu.help    = {''
'Choose unwarped image volumes acquired with blip-up.'
};
hysco_bu.filter  = 'image';
hysco_bu.ufilter = '.*';
hysco_bu.num     = [0 Inf];

% unwarped blip-dw images
hysco_bd         = cfg_files;
hysco_bd.tag     = 'hysco_bd';
hysco_bd.name    = 'Unwarped blip-down images';
hysco_bd.help    = {''
'Choose unwarped image volumes acquired with blip-down.'
};
hysco_bd.filter  = 'image';
hysco_bd.ufilter = '.*';
hysco_bd.num     = [0 Inf];


% % define weights for blip-up images
% weights_bu      = others_bu;
% weights_bu.tag  = 'weights_bu';
% weights_bu.name = 'Weights for blip-up images';
% weights_bu.val  = {{''}};
% weights_bu.help = {
%     'Weights for weighted combination.'};

% % define weights for blip-down images
% weights_bd      = others_bd;
% weights_bd.tag  = 'weights_bd';
% weights_bd.name = 'Weights for blip-down images';
% weights_bd.val  = {{''}};
% weights_bd.help = {
%     'Weights for weighted combination.'};

% multihysco_write apply susceptibility correction
dummy_wcomb        = cfg_menu;
dummy_wcomb.tag    = 'dummy_wcomb';
dummy_wcomb.name   = 'Write combined images';
dummy_wcomb.help   = {''
'Specify what is written out.'
''};
dummy_wcomb.labels = {'Write nothing', 'Write weighted average', 'Write arithmetic average', 'Write both weighted and arithmetic average'};
dummy_wcomb.values = {0 1 2 3};
dummy_wcomb.val    = {1};

% is the fieldmap from FSL or ACID?
dummy_fmaptype        = cfg_menu;
dummy_fmaptype.tag    = 'dummy_fmaptype';
dummy_fmaptype.name   = 'ACID or FSL fieldmap';
dummy_fmaptype.help   = {'Specify the toolbox from which the fieldmap was generated. By default, an ACID fieldmap is assumed (the name starts with "HySCO").'};
dummy_fmaptype.labels = {'ACID','FSL'};
dummy_fmaptype.values = {1 0};
dummy_fmaptype.val    = {1};

% Exbranch
hysco_wcomb_multi      = cfg_exbranch;
hysco_wcomb_multi.tag  = 'hysco_wcomb_multi';
hysco_wcomb_multi.name = 'HySCO: combine blip-up and blip-down images';
hysco_wcomb_multi.val  = {source_bu_apply hysco_bu hysco_bd fmap dummy_fmaptype dummy_wcomb phdir};
hysco_wcomb_multi.help = {
                    'This option applies the field inhomogeneity maps from the Hyperelastic susceptibility artifact correction of diffusion weighted images to undistort images of choice.'
                    'The output is either the arithmetic mean or the weighted combination.'
};
hysco_wcomb_multi.prog = @local_hysco_comb;
hysco_wcomb_multi.vout = @vout_hysco_comb;

%% HYSCO choice
hysco_choice         = cfg_choice;
hysco_choice.tag     = 'hysco_choice';
hysco_choice.name    = 'HySCO';
hysco_choice.help    = {
                 'Estimate field-inhomogeneities from blip-up and -down images and apply it to undistort the images.'
                 'Apply field-inhomogneities to undistorted other image(s).'
}';
if dummy_hysco_apply_multi
    hysco_choice.values = {hysco hysco_apply hysco_wcomb_multi};
else
    hysco_choice.values = {hysco hysco_apply};
end

end

function out = local_hysco(job)

    alpha       = acid_get_defaults('hysco.alpha');
    beta        = acid_get_defaults('hysco.beta');
    restrictdim = acid_get_defaults('hysco.restrictdim');
    dummy_ecc   = acid_get_defaults('hysco.dummy_ecc');
    res         = acid_get_defaults('hysco.resample');

    % run HySCO
    [NiB, Ni1, Ni2] = acid_hysco(char(job.source_bu), char(job.source_bd), char(job.others_bu), char(job.others_bd), ...
        job.phdir, job.dummy_fast, dummy_ecc, alpha, beta, restrictdim, res);

    % get field inhomogeneity map
    out.fieldmap = {NiB.dat.fname};

    % get saved HySCO corrected images - blip up
    if ~isempty(Ni1)
        out.bu_others = {[Ni1.dat.fname ',1']};
    end
    
    % get saved HySCO corrected images - blip down
    if ~isempty(Ni2)
        out.bd_others = {[Ni2.dat.fname ',2']};
    end
end

function dep = vout_hysco(~)
    kk = 1;
    dep(1)             = cfg_dep;
    dep(1).sname       = 'Unwarped blip-up images';
    dep(1).src_output  = substruct('.','bu_others');
    dep(1).tgt_spec    = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk+1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Unwarped blip-down images';
    dep(kk).src_output = substruct('.','bd_others');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk + 1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Field inhomogeneity map';
    dep(kk).src_output = substruct('.','fieldmap');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

function out = local_hysco_apply(job)
    [V_bu, V_bd] = acid_hysco_apply(char(job.source_bu), char(job.others_bu), char(job.others_bd), char(job.fmap), job.phdir);

    % get field inhomogeneity map
    out.fieldmap = job.fmap(:);
    
    % get saved HySCO corrected images - blip up
    if logical(~isempty(V_bu))
        out.bu_others = {[V_bu.fname ',1']};
    elseif logical(isempty(V_bu))
        out.bu_others = [];
    end

    % get saved HySCO corrected images - blip down
    if logical(~isempty(V_bd))
        out.bd_others = {[V_bd.fname ',1']};
    elseif logical(isempty(V_bd))
        out.bd_others = [];
    end  

end

function dep = vout_hysco_apply(~)
    kk=1;
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Unwarped blip-up images';
    dep(1).src_output = substruct('.','bu_others');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk+1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Unwarped blip-down images';
    dep(kk).src_output = substruct('.','bd_others');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk + 1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Field inhomogeneity map';
    dep(kk).src_output = substruct('.','fieldmap');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

function out = local_hysco_write_multi(job)
    [~,V14D,V24D,VJacup,VJacdw] = acid_hysco_write_multiB(char(job.source_bu), char(job.others_up), char(job.others_dw), char(job.fmap), job.phdir, 1);

    out.fieldmap   = job.fmap(:);
    if logical(~isempty(VJacup))
        out.rothers_jacup = {VJacup(1).fname};
    else
        out.rothers_jacup = [];
    end
    if logical(~isempty(VJacdw))
        out.rothers_jacdw = {VJacdw(1).fname};
    else
        out.rothers_jacdw = [];
    end

    if logical(~isempty(V14D))
        out.rothers_up = {V14D(1).fname};
    elseif(logical(isempty(V14D)))
        out.rothers_up = [];
    end
    out.rothers_up = out.rothers_up(:);

    if logical(~isempty(V24D))
        out.rothers_dw = {V24D(1).fname};
    elseif(logical(isempty(V24D)))
        out.rothers_dw = [];
    end
    out.rothers_dw = out.rothers_dw(:);

end

function dep = vout_hysco_write_multi(~)
    kk=1;
    dep(1)            = cfg_dep;
    dep(1).sname      = 'Unwarped blip-up images';
    dep(1).src_output = substruct('.','rothers_up');
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk+1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Unwarped blip-down images';
    dep(kk).src_output = substruct('.','rothers_dw');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk+1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Field inhomogeneity map';
    dep(kk).src_output = substruct('.','fieldmap');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk + 1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Estimated Jacoby up';
    dep(kk).src_output = substruct('.','rothers_jacup');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    kk = kk + 1;
    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Estimated Jacoby dw';
    dep(kk).src_output = substruct('.','rothers_jacdw');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

function out = local_hysco_comb(job)
        
    [Ni_aa, Ni_wa] = acid_hysco_comb(char(job.source_bu), char(job.hysco_bu), char(job.hysco_bd), char(job.fmap), job.dummy_fmaptype, job.dummy_wcomb, job.phdir);

    out.fieldmap = job.fmap(:);
    
    if job.dummy_wcomb == 2 || job.dummy_wcomb == 3
        out.comb_aa{1,:} = [Ni_aa(1).dat.fname ',1'];
    end

    if job.dummy_wcomb == 1 || job.dummy_wcomb == 3
        out.comb_wa{1,:} = [Ni_wa(1).dat.fname ',1'];
    end
end

function dep = vout_hysco_comb(job)

    kk = 0;
    if job.dummy_wcomb == 1 || job.dummy_wcomb == 3
        kk = kk + 1;
        dep(1)            = cfg_dep;
        dep(1).sname      = 'Weighted average of blip-up and blip-down images';
        dep(1).src_output = substruct('.','comb_wa');
        dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        kk = kk + 1;
    end

    if job.dummy_wcomb == 2 || job.dummy_wcomb == 3
        if job.dummy_wcomb == 2
            kk = kk + 1;   
        end
        dep(kk)            = cfg_dep;
        dep(kk).sname      = 'Arithmetic average of blip-up and blip-down images';
        dep(kk).src_output = substruct('.','comb_aa');
        dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        kk = kk + 1;
    end

    if job.dummy_wcomb == 0
        kk = kk + 1;   
    end

    dep(kk)            = cfg_dep;
    dep(kk).sname      = 'Estimated fieldmap';
    dep(kk).src_output = substruct('.','fieldmap');
    dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
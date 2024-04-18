function [V14D, V24D] = acid_hysco_apply(P_ref, P_bu, P_bd, P_field, phdir)

% =========================================================================
% (c) Lars Ruthotto and Siawoosh Mohammadi 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
% 
% FORMAT acid_hysco_write(POI1,POI2,Bc,pe_direction)
%
% Applying HySCO result to other blip-up or  driver for HySCO (Hyperelastic Susceptibility COrrection of DTI)
%
% Inputs:
%   P_ref   - reference blip-up volume
%   P_bu    - filenames of additional blip-up volumes
%   P_bd    - filenames of additional blip-down volumes
%   P_field - filename of inhomogeneity estimate produced by HySCO
%   phdir   - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
%                 (data dimensions will be flipped accordingly)
%   p_out   - output directory
%
% Please cite one of the following works when using this software.
%
% @inproceedings{Ruthotto2013,
%   author    = {Ruthotto, L and Mohammadi, S and Heck, C and Modersitzki, J and Weiskopf, N},
%   title     = {HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM}},
%   booktitle = {Bildverarbeitung f{\"u}r die Medizin 2013},
%   year      = {2013}
% }
%
% @article{Ruthotto2012,
%   author  = {Ruthotto, L and Kugel, H and Olesch, J and Fischer, B and Modersitzki, J and Burger, M and Wolters, CH},
%   title   = {Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images}},
%   journal = {Physics in Medicine and Biology},
%   volume  = {57},
%   number  = {18},
%   pages   = {5715--5731}
%   year    = {2012}
% }
%
    %{
        (c) Lars Ruthotto and Jan Modersitzki 2013

        This file is part of HySCO (Version 1.0, 2013/03/28)
                               -  Hyperelastic Susceptibility Artefact Correction for DTI


        HySCO is free but copyright software, distributed under the terms of the 
        GNU General Public Licence as published by the Free Software Foundation 
        (Version 3, 29 June 2007) http://www.gnu.org/licenses/gpl.html


        This code is provided "as is", without any warranty of any kind, either
        expressed or implied, including but not limited to, any implied warranty
        of merchantibility or fitness for any purpose. In no event will any party
        who distributed the code be liable for damages or for any claim(s) by
        any other party, including but not limited to, any lost profits, lost
        monies, lost data or data rendered inaccurate, losses sustained by
        third parties, or any other special, incidental or consequential damages
        arising out of the use or inability to use the program, even if the
        possibility of such damages has been advised against. The entire risk
        as to the quality, the performace, and the fitness of the program for any
        particular purpose lies with the party using the code.

        This code is especially not intended for any clinical or diagnostic use. 
    %}

% =========================================================================

VG   = spm_vol(P_ref);
V_bu = spm_vol(P_bu);
V_bd = spm_vol(P_bd);

if ~exist('res','var')
    res = -4;
end

V14D = [];
V24D = [];

% extract data resolution and domain info 
% Note that domain is assumed to be rectangular and aligned to the coordinate system,
% i.e. omega = [omega(1),omega(2)] x [omega(3),omega(4)] x [omega(5),omega(6)]
if numel(V_bu)>=1
    m    = VG.dim;
    Vmat = sqrt(sum(VG.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
elseif numel(V_bd)>=1
    m    = VG.dim;
    Vmat = sqrt(sum(VG.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
else
    error('No input volumes supplied!');
end
    
omega = zeros(1,6);
omega(2:2:end) = Vmat(1:3).*m; % modified by SM to make sure that voxel size is kept

% load in 4D images
V_bu = acid_load_4Dimage(P_bu);
V_bd = acid_load_4Dimage(P_bd);

% define output directory
keyword = 'HySCO';
[path,fname,~] = spm_fileparts(VG.fname);

p_out = acid_bids(path,fname,keyword,1);


% find out if inhomogeneity is nodal or staggered
Bc      = spm_read_vols(spm_vol(P_field)); 
isNodal = numel(Bc)==prod(m+1);
mstg    = m; mstg(phdir) = mstg(phdir)+1;
isStg   = numel(Bc)==prod(mstg);
if not(isNodal) && not(isStg)
    error('resolution of inhomogeneity not compatible with data');
end

% permute data dimensions such that phase encoding is along second index
if isNodal
    switch phdir
        case 1
            read_data  = @(str) permute(sm_read_vols(str,VG,res),[2 1 3]);
            write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(ipermute(A,[2 1 3]), V, VG, Ni, vol, n_vols, dtype, fname_out);
            omega = omega([3 4 1 2 5 6]);  m = m([2 1 3]);
            vecperm = [2 1 3];
        case 2
            read_data  = @(str) sm_read_vols(str,VG,res);
            write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(A, V, VG, Ni, vol, n_vols, dtype, fname_out);
            vecperm = [1 2 3];
        case 3
            read_data  = @(str) permute(sm_read_vols(str,VG,res),[1 3 2]);
            write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(ipermute(A,[1 3 2]), V, VG, Ni, vol, n_vols, dtype, fname_out);
            omega = omega([1 2 5 6 3 4]);  m = m([1 3 2]);
            vecperm = [1 3 2];
    end
elseif isStg
    switch phdir
        case 1
            read_data  = @(str) sm_read_vols(str,VG(1),res);
            write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(A, V, VG, Ni, vol, n_vols, dtype, fname_out);
            vecperm = [1 2 3];
        case 2
            read_data  = @(str) permute(sm_read_vols(str,VG(1),res),[2 1 3]);
            write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(ipermute(A,[2 1 3]), V, VG, Ni, vol, n_vols, dtype, fname_out);
            omega = omega([3 4 1 2 5 6]);  m = m([2 1 3]);
            vecperm = [2 1 3];
        case 3
            read_data  = @(str) permute(sm_read_vols(str,VG(1),res),[3 1 2]);
            write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(ipermute(A,[3 1 2]), V, VG, Ni, vol, n_vols, dtype, fname_out);
            omega = omega([ 5 6 1:4]);  m = m([3 1 2]);
            vecperm = [3 1 2];
    end
end

dim = numel(find(m>1));
if dim==2
    if m(3)>1
        error('Invalid pe direction for 2D image. Number of voxels is m=[%d,%d]',m);
    end
    m = m(1:2);
    omega= omega(1:4);
end

Bc = permute(spm_read_vols(spm_vol(P_field)),vecperm);
h  = (omega(2:2:end)-omega(1:2:end))./m;

% compute transformations and intensity modulations
if isNodal
    y1  = acid_hysco_getTrafoEPI(Bc,[0;1;0],omega,m,'matrixFree',1);
    y2  = acid_hysco_getTrafoEPI(Bc,[0;-1;0],omega,m,'matrixFree',1);
    y1  = nodal2center(y1,m);
    y2  = nodal2center(y2,m);
    pB  = acid_hysco_getPartialB(Bc,omega,m,'cc','matrixFree',1);
elseif isStg
    xc  = reshape(getCellCenteredGrid(omega,m),[],3);
    Bc  = reshape(Bc,m+eye(1,dim));               % 1-staggered
    Bcc = .5*(Bc(1:end-1,:,:) + Bc(2:end,:,:));   % cell-centered
    pB  = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
    y1  = xc; y1(:,1) = y1(:,1) + Bcc(:);
    y2  = xc; y2(:,1) = y2(:,1) - Bcc(:);  
end
Jac1 = 1 + pB;
Jac2 = 1 - pB;

fname_bd = acid_bids_filename(V_bd(1), keyword, '_dwi', '.nii');
fname_bu = acid_bids_filename(V_bu(1), keyword, '_dwi', '.nii');
fname_bd = [p_out filesep fname_bd];
fname_bu = [p_out filesep fname_bu];

% apply field inhomogeneity map
for vol = 1:numel(V_bu)
    fprintf('Apply estimated field inhomogeneity to blip-up volumes %d\n',vol);
   % I1 = read_data(V_bu(vol));
   % I1opt = reshape(linearInterMex(I1,omega,y1).*Jac1(:),m);
    I1    = getSplineCoefficients(read_data(V_bu(vol)));
    I1opt = reshape(splineInterMex(I1,omega,y1).*Jac1(:),m);
    
    if vol==1
        Ni1 = write_data(I1opt, V_bu(vol), VG, [], vol, numel(V_bu), 'uint16', fname_bu);
    else
        Ni1 = write_data(I1opt, V_bu(vol), VG, Ni1, vol, numel(V_bu), 'uint16', fname_bu);
    end
end

for vol = 1:numel(V_bd)
    fprintf('Apply estimated field inhomogeneity to blip-down volumes %d\n',vol);
   % I2 = read_data(V_bd(vol));
   % I2opt = reshape(linearInterMex(I2,omega,y2).*Jac2(:),m);
    I2    = getSplineCoefficients(read_data(V_bd(vol)));
    I2opt = reshape(splineInterMex(I2,omega,y2).*Jac2(:),m);
    
    if vol==1
        Ni2 = write_data(I2opt, V_bd(vol), VG, [], vol, numel(V_bd), 'uint16', fname_bd);
    else
        Ni2 = write_data(I2opt, V_bd(vol), VG, Ni2, vol, numel(V_bd), 'uint16', fname_bd);
    end
end

% save json files
if ~isempty(V_bu), acid_save_json(V_bu(1), p_out, keyword); end
if ~isempty(V_bd), acid_save_json(V_bd(1), p_out, keyword); end
    
% save bvals and bvecs files
if ~isempty(V_bu), acid_save_bvals_bvecs(V_bu(1), p_out, keyword); end
if ~isempty(V_bd), acid_save_bvals_bvecs(V_bd(1), p_out, keyword); end
    
end

function Atmp = sm_read_vols(strS,strT,res)
    V    = strS;
    VG   = strT;
    Atmp = acid_read_vols(V,VG,res);
end

function Ni = spm_write_image_4d(I, V, VG, Ni, vol, n_vols, dtype, fname_out)
    
    if ~exist('VG','var')
        VG = V;
    end
    if isempty(VG)
        VG = V;
    end

    if vol == 1 
        Ni      = nifti;
        Ni.mat  = VG(1).mat;
        Ni.mat0 = VG(1).mat;
        
        if n_vols > 1
            dm = [VG(1).dim n_vols];
            Ni.descrip = '4d array of HYSCO corrected images';
        else
            dm = VG(1).dim;
        end
        
        Ni.dat = file_array(fname_out, dm, dtype, 0, 1, 0);
        create(Ni);
        
        if n_vols > 1
            spm_progress_bar('Init', n_vols, Ni.descrip, 'volumes completed');
        end      
    end

    % select the image to write
    Ni.dat(:,:,:,vol) = I;
    spm_progress_bar('Set',vol);

    disp(['Image #: ' num2str(vol) ' undistorted'])
    if vol == size(V,1)
        spm_progress_bar('Clear');
    end
end
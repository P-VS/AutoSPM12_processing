function [Ni_aa, Ni_wa] = acid_hysco_comb(P_ref, P_bu, P_bd, P_fmap, dummy_fmaptype, dummy_wcomb, phdir)

% =========================================================================
% Inputs:
%  P_ref        - filename of reference blip-up volume
%  P_bu         - filenames for additional blip-up volumes
%  P_bd         - filenames for additional blip-down volumes
%  P_weights_bu - filenames of weights for blip-up volumes (optional)
%  P_weights_bd - filenames of weights for blip-down volumes (optional)
%  P_fname      - filename of inhomogeneity estimate produced by HySCO
%  phdir        - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
%                 (data dimensions will be flipped accordingly)
%  dummy_3dor4d - 3D or 4D dataset
%  dummy_wcomb  - write data (arithmetic mean, weighted combination, or both)
%  dummy_fmaptype   - the fieldmap is based on HySCO or FSL topup
%
% Please cite one of the following works when using this software
%
% Based on acid_hysco to allow for combination of multiple HySCO
% fieldmaps as well as their weighted combination.
% S. Mohammadi 2.10.2019
% =========================================================================
tic
VG   = spm_vol(P_ref);
V_bu = spm_vol(P_bu);
V_bd = spm_vol(P_bd);

% default parameters
if ~exist('res','var')
    res = -7;
end

if ~exist('kt','var')
    if dummy_fmaptype==1
        kt = 10;
    else
        kt = 5;%20;
    end
end

if ~exist('k0','var')
    k0 = 1;
end

Vw1 = [];
Vw2 = [];
Ni_aa = [];
Ni_wa = [];

% extract data resolution and domain info
% Note that domain is assumed to be rectangular and aligned to the coordinate system,
% i.e. omega = [omega(1),omega(2)] x [omega(3),omega(4)] x [omega(5),omega(6)]
if numel(V_bu)>=1
    m     = VG.dim;
    Vmat  = sqrt(sum(VG.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
elseif numel(V_bd)>=1
    m     = VG.dim;
    Vmat  = sqrt(sum(VG.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
else
    error('No input volumes supplied!');
end

omega = zeros(1,6);
omega(2:2:end) = Vmat(1:3).*m; % modified by SM to make sure that voxel size is kept

% load in 4D images
V_bu = acid_load_4Dimage(P_bu);
V_bd = acid_load_4Dimage(P_bd);

% define output directory
keyword = 'HySCO-COMB';
[path,fname,~] = spm_fileparts(VG.fname);

p_out = acid_bids(path,fname,keyword,1);


diary([p_out filesep 'logfile_' keyword '.txt'])

% permute data dimensions such that phase encoding is along first index
if phdir == 1
    fprintf('\nPhase-encoding direction: x\n')
elseif phdir == 2
    fprintf('\nPhase-encoding direction: y\n')
elseif phdir == 3
    fprintf('\nPhase-encoding direction: z\n')
end


if ~isempty(P_fmap)
    
    % find out if inhomogeneity is nodal or staggered
    Bc = spm_read_vols(spm_vol(P_fmap(1,:)));
    
    if dummy_fmaptype % ACID field map
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
                    write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,Ni,Nvol,sz,dtype,fname_out);
                    omega = omega([3 4 1 2 5 6]);  m= m([2 1 3]);
                    vecperm = [2 1 3];
                case 2
                    read_data   = @(str) sm_read_vols(str,VG,res);
                    write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(A,V,VG,Ni,Nvol,sz,dtype,fname_out);
                    vecperm = [1 2 3];
                case 3
                    read_data   = @(str) permute(sm_read_vols(str,VG,res),[1 3 2]);
                    write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(ipermute(A,[1 3 2]),V,VG,Ni,Nvol,sz,dtype,fname_out);
                    omega = omega([1 2 5 6 3 4]);  m=m([1 3 2]);
                    vecperm = [1 3 2];
            end
        elseif isStg
            switch phdir
                case 1
                    read_data   = @(str) sm_read_vols(str,VG(1),res);
                    write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(A,V,VG,Ni,Nvol,sz,dtype,fname_out);
                    vecperm = [1 2 3];
                case 2
                    read_data   = @(str) permute(sm_read_vols(str,VG(1),res),[2 1 3]);
                    write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,Ni,Nvol,sz,dtype,fname_out);
                    omega = omega([3 4 1 2 5 6]);  m= m([2 1 3]);
                    vecperm = [2 1 3];
                case 3
                    read_data   = @(str) permute(sm_read_vols(str,VG(1),res),[3 1 2]);
                    write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(ipermute(A,[3 1 2]),V,VG,Ni,Nvol,sz,dtype,fname_out);
                    omega = omega([ 5 6 1:4]);  m=m([3 1 2]);
                    vecperm = [3 1 2];
            end
        end
        
        % compute inhomogeneity
        Bc = zeros(size(permute(spm_read_vols(spm_vol(P_fmap(1,:))),vecperm)));
        for inx = 1:size(P_fmap,1)
            Bctmp = permute(spm_read_vols(spm_vol(P_fmap(inx,:))),vecperm);
            Bc = Bctmp + Bc;
        end

        % compute transformations and intensity modulations
        if isNodal
            y1    = acid_hysco_getTrafoEPI(Bc,[0;1;0],omega,m,'matrixFree',1);
            y2    = acid_hysco_getTrafoEPI(Bc,[0;-1;0],omega,m,'matrixFree',1);
            y1    = nodal2center(y1,m);
            y2    = nodal2center(y2,m);
            pB    = acid_hysco_getPartialB(Bc,omega,m,'cc','matrixFree',1);
        elseif isStg
            xc    = reshape(getCellCenteredGrid(omega,m),[],3);
            Bc    = reshape(Bc,m+[1,0,0]);                  % 1-staggered
            Bcc   = .5*(Bc(1:end-1,:,:) + Bc(2:end,:,:));   % cell-centered
            y1    = xc; y1(:,1) = y1(:,1) + Bcc(:);
            y2    = xc; y2(:,1) = y2(:,1) - Bcc(:);
            h     = (omega(2:2:end)-omega(1:2:end))./m;
            pB    = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
        end
    elseif dummy_fmaptype==0 % FSL field map
        AB = acid_read_vols(spm_vol(P_fmap),VG,1);
        switch phdir
            case 1
                read_data   = @(str) permute(sm_read_vols(str,VG,res),[2 1 3]);
                write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[2 1 3]),V,prefix,VG,dummy_bids,keyword,path_hysco);
                write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,prefix,Ni,Nvol,sz,dummy_bids,keyword,path_hysco);
                omega=omega([3 4 1 2 5 6]);  m= m([2 1 3]);
                vecperm = [2 1 3];
            case 2
                read_data   = @(str) sm_read_vols(str,VG,res);
                write_data  = @(A,V,prefix) spm_write_image(A,V,prefix,VG,dummy_bids,keyword,path_hysco);
                write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(A,V,VG,prefix,Ni,Nvol,sz,dummy_bids,keyword,path_hysco);
                vecperm = [1 2 3];
            case 3
                read_data   = @(str) permute(sm_read_vols(str,VG,res),[1 3 2]);
                write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[1 3 2]),V,prefix,VG,dummy_bids,keyword,path_hysco);
                write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[1 3 2]),V,VG,prefix,Ni,Nvol,sz,dummy_bids,keyword,path_hysco);
                omega=omega([1 2 5 6 3 4]);  m=m([1 3 2]);
                vecperm = [1 3 2];
        end
        h             = (omega(2:2:end)-omega(1:2:end))./m;
        Bc            = permute(AB,vecperm);
        pB            =  permute(zeros(m),vecperm);
        pB(2:end,:,:) = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
    end
    
    % compute Jacobians
    Jac1 = 1 + pB; Jac1(Jac1<0)=0;
    Jac2 = 1 - pB; Jac2(Jac2<0)=0;
    
    % save Jacobians
    fname = acid_bids_filename(VG, ['HySCO' '-Jac1'], '', '.nii');
    fname_out = [p_out filesep fname]; 
    write_data(Jac1, VG, [], [], 1, 1, 'float32', fname_out);
    
    fname = acid_bids_filename(VG, ['HySCO' '-Jac2'], '', '.nii');
    fname_out = [p_out filesep fname]; 
    write_data(Jac2, VG, [], [], 1, 1, 'float32', fname_out);

    % write weighted average of blip-up and blip-down images
    if dummy_wcomb==1 || dummy_wcomb==3
             
        for vol = 1:numel(V_bu)
            fprintf('Estimate weighted combination for volume %d\n',vol);
            I1 = read_data(V_bu(vol));
            I2 = read_data(V_bd(vol));
            
            if dummy_fmaptype
                f1 = 2-acid_hysco_fermi(Jac1,k0,kt); % we take 2- to downweigh regions that were squeezed
                f2 = 2-acid_hysco_fermi(Jac2,k0,kt); % we take 2- to downweigh regions that were squeezed
            else
                f1 = acid_hysco_fermi(Jac1,k0,kt); % we take acid_hysco_fermi to downweigh regions that were squeezed
                f2 = acid_hysco_fermi(Jac2,k0,kt); % we take acid_hysco_fermi to downweigh regions that were squeezed
            end
            
            if ~logical(isempty(Vw1)) && ~logical(isempty(Vw2))
                I_w1 = read_data(Vw1(vol));
                I_w1(I_w1==0)=1;
                I_w2 = read_data(Vw2(vol));
                I_w2(I_w2==0)=1;
            else
                I_w1 = 1;
                I_w2 = 1;
            end
            
            I_wa = (I1.*f1.*I_w1 + I2.*f2.*I_w2)./(f1.*I_w1+f2.*I_w2);
            fprintf('Write weighted average of blip-up and blip-down volumes %d\n',vol);
            
            % save combined data
            fname = acid_bids_filename(V_bu(vol), ['COMB' '-WA'], '_dwi', '.nii');
            fname_out = [p_out filesep fname];
            
            if vol==1
                Ni_wa = write_data(I_wa, V_bu(vol), VG, [], vol, numel(V_bu), 'uint16', fname_out);
            else
                Ni_wa = write_data(I_wa, V_bu(vol), VG, Ni_wa, vol, numel(V_bu), 'uint16', fname_out);
            end
        end
    end
end

% write arithmetic average of blip-up and blip-down images
if dummy_wcomb==2 || dummy_wcomb==3
    if isempty(P_fmap)
        switch phdir
            case 1
                read_data  = @(str) permute(sm_read_vols(str,VG,res),[2 1 3]);
                write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,Ni,Nvol,sz,dtype,fname_out);
            case 2
                read_data  = @(str) sm_read_vols(str,VG,res);
                write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(A,V,VG,Ni,Nvol,sz,dtype,fname_out);
            case 3
                read_data  = @(str) permute(sm_read_vols(str,VG,res),[1 3 2]);
                write_data = @(A,V,VG,Ni,Nvol,sz,dtype,fname_out) spm_write_image_4d(ipermute(A,[1 3 2]),V,VG,Ni,Nvol,sz,dtype,fname_out);
        end
    end

    % write arithmetic average of blipup and down image volumes
    for vol = 1:numel(V_bu)
        fprintf('Estimate weighted combination for volume %d\n',vol);
        I1 = read_data(V_bu(vol));
        I2 = read_data(V_bd(vol));

        I_aa = (I1 + I2)*0.5;
        
        % save combined data
        fname = acid_bids_filename(V_bu(vol), ['COMB' '-AA'], '_dwi', '.nii');
        fname_out = [p_out filesep fname];
        
        if vol==1
            Ni_aa = write_data(I_aa, V_bu(vol), VG, [], vol, numel(V_bu), 'uint16', fname_out);
        else
            Ni_aa = write_data(I_aa, V_bu(vol), VG, Ni_aa, vol, numel(V_bu), 'uint16', fname_out);
        end  
    end
end

T = toc/60;

T = duration(minutes(T),'format','hh:mm:ss');
disp(['The total time for ' keyword ' was: ' char(T) '.']);
diary off


end

function Ni = spm_write_image_4d(I, V, VG, Ni, Nvol, sz, dtype, fname_out)
    
    if ~exist('VG','var')
        VG = V;
    end
    
    if isempty(VG)
        VG = V;
    end

    if Nvol == 1 
        Ni      = nifti;
        Ni.mat  = VG(1).mat;
        Ni.mat0 = VG(1).mat;
        
        if sz > 1
            dm = [VG(1).dim sz];
            Ni.descrip = '4d array of HYSCO corrected images';
        else
            dm = VG(1).dim;
        end
        
        Ni.dat = file_array(fname_out, dm, dtype, 0, 1, 0);
        create(Ni);
        
        if sz > 1
            spm_progress_bar('Init', sz, Ni.descrip, 'volumes completed');
        end      
    end

    % select the image to write
    Ni.dat(:,:,:,Nvol) = I;
    spm_progress_bar('Set',Nvol);

    disp(['Image #: ' num2str(Nvol) ' undistorted'])
    if Nvol == size(V,1)
        spm_progress_bar('Clear');
    end
end


function Atmp = sm_read_vols(strS,strT,res)
% =========================================================================
% function sm_write_image(strS,strT,res)
%
% read image volume
%
% Input:
%   str1   - structure containing image volume information of ith image
%   V      - structure containing image volume information of target image
%   res    - resampling order
%   perm   - permutation
%
% =========================================================================

    V       = strS;
    VG      = strT;
    Atmp    = acid_read_vols(V,VG,res);
    Atmp(isnan(Atmp(:)))=0;
    
end
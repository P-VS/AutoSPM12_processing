function [NiB, Ni1, Ni2] = acid_hysco(P_bu_ref, P_bd_ref, P_bu, P_bd, phdir, fullres, doECC, alpha, beta, restrictdim, res)

% =========================================================================
% (c) Lars Ruthotto and Siawoosh Mohammadi 2016
% http://www.mathcs.emory.edu/~lruthot/
%
% function acid_hysco_main(PI1,PI2,POI1,POI2,pe_direction,full_res,doECC,alpha,beta)
%
% Main driver for HySCO 2.0 (Hyperelastic Susceptibility COrrection of DTI)
%
% The inexact Gauss-Newton method with block Jacobi preconditioner as
% described in~\cite{MacdonaldRuthotto2016} is used.
%
% We thank Jan Macdonald for his help developing the new optimization.
%
% Inputs:
%
%  P_bu_ref     - filename of reference blip-up
%  P_bd_ref     - filename of reference blip-down
%  P_bu         - matrix of filenames for additional blip-up volumes
%  P_bd         - matrix of filenames for additional blip-down volumes
%  phdir        - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
%                 (data dimensions will be flipped accordingly)
%  full_res     - finest level for multi-level correction, boolean
%  doECC        - do nonlinear eddy-correction for additinal volumes (equal
%                 number of blip-up and blip-down data required.)
%  alpha        - regularization paramter for diffusion term
%  beta         - regularization paramter for Jacobian term
%
% Please cite one of the following works when using this software
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
% @article{MacdonaldRuthotto2016,
%   author        = {Macdonald, J and Ruthotto, L},
%   title         = {Efficient Numerical Optimization For Susceptibility Artifact Correction Of EPI-MRI}},
%   archivePrefix = {arXiv},
%   eprint        = {xxxx.xxxx},
%   primaryClass  = {xx-xx},
%   year          = {2016}
% }
%
%{
    (c) Lars Ruthotto 2016

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

tic

plots = 0;
doReportResults = 0;

Ni1 = [];
Ni2 = [];

if isempty(restrictdim)
    restrictdim = [1,1,1];
end

V_bu_ref = spm_vol(P_bu_ref);
if size(V_bu_ref,1)>1
    V_bu_ref = V_bu_ref(1);
end

V_bd_ref = spm_vol(P_bd_ref);
if size(V_bd_ref,1)>1
    V_bd_ref = V_bd_ref(1);
end

% load in 4D images

V_bu = acid_load_4Dimage(P_bu);
V_bd = acid_load_4Dimage(P_bd);

if isempty(V_bu) && ~isempty(V_bd)
    up_or_down_struct = V_bd(1,:);
elseif ~isempty(V_bu)
    up_or_down_struct = V_bu(1,:);
else
    up_or_down_struct = V_bu_ref;
end

% define output directory
keyword = 'HySCO';
[path,fname,~] = spm_fileparts(up_or_down_struct.fname);

p_out = acid_bids(path,fname,keyword,1);

% start logging
diary([p_out filesep 'logfile_' keyword '.txt'])

% automatically generate HTML report of correction results (not tested yet)
if doReportResults
    Report = @(varargin) acid_hysco_report(varargin{:});
else
    Report = @(varargin) [];
end

% write header of report
reportDir  = [p_out filesep 'reports'];
reportName = fname;
reportFile = Report('header',reportDir,reportName);

% extract data resolution and domain info
% Note that domain is assumed to be rectangular and aligned to the coordinate system,
% i.e. omega = [omega(1),omega(2)] x [omega(3),omega(4)] x [omega(5),omega(6)]

m = V_bu_ref(1).dim;
% check if all blip-up and blip-down images are of same resolution
if any(V_bd_ref(1).dim~=m); error('%s: blip-up and blip-down images must have same resolution',mfilename); end
for vol = 1:numel(V_bu)
    if any(V_bu(vol).dim ~=m)
        error('%s: dimensions of reference blip-up and other blip-up images must match! Violated at least for vol=%d!',mfilename,vol)
    end
end
for vol = 1:numel(V_bd)
    if any(V_bd(vol).dim ~=m)
        error('%s: dimensions of reference blip-down and other blip-down images must match! Violated at least for vol=%d!',mfilename,vol)
    end
end

omega   = zeros(1,6);
Vmat    = sqrt(sum(V_bu_ref(1).mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
omega(2:2:end) = Vmat(1:3).*m; % modified by SM to make sure that voxel size is kept
h       = (omega(2:2:end)-omega(1:2:end))./m;
VG      = V_bu_ref(1);

% permute data dimensions such that phase encoding is along first index
if phdir == 1
    fprintf('\nPhase-encoding direction: x\n')
elseif phdir == 2
    fprintf('\nPhase-encoding direction: y\n')
elseif phdir == 3
    fprintf('\nPhase-encoding direction: z\n')
end



switch phdir
    case 1
        read_data  = @(str) sm_read_vols(str,V_bu_ref(1),res);
        write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(A, V, VG, Ni, vol, n_vols, dtype, fname_out);
    case 2
        alpha      = permute(alpha,[2 1 3]);
        read_data  = @(str) permute(sm_read_vols(str,V_bu_ref(1),res),[2 1 3]);
        write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(ipermute(A,[2 1 3]), V, VG, Ni, vol, n_vols, dtype, fname_out);
        omega = omega([3 4 1 2 5 6]);
        m = m([2 1 3]);
    case 3
        alpha      = permute(alpha,[3 1 2]);
        read_data  = @(str) permute(sm_read_vols(str,V_bu_ref(1),res),[3 1 2]);
        write_data = @(A, V, VG, Ni, vol, n_vols, dtype, fname_out) spm_write_image_4d(ipermute(A,[3 1 2]), V, VG, Ni, vol, n_vols, dtype, fname_out);
        omega=omega([ 5 6 1:4]);
        m = m([3 1 2]);
end

dim = numel(find(m>1));
if dim==2
    if m(3)>1
        error('Invalid pe direction for 2D image. Number of voxels is m=[%d,%d]',m);
    end
    m = m(1:2);
    omega= omega(1:4);
end

% normalize intensities for field estimation
[I1, I2] = normalizeIntensities(read_data(V_bu_ref(1)),read_data(V_bd_ref(1)));

% get a suitable discretization size for multi-level, which is divisible by 4 and
% such that on the coarsest level has at least 4 cells in each direction
[MLdata,minLevel,maxLevel] = getMultilevel({I1,I2},omega,m,'fig',0,'restrictdim',restrictdim);
minLevel = max(maxLevel-2,minLevel);
if not(fullres) % perform correction only until half resolution (if possible)
    maxLevel = max(minLevel,maxLevel - 1);
end

% configure FAIR's viewer module
viewImage('reset','viewImage','imgmontage','colormap','gray');

% configure cubic B-spline interpolation with moments regularization
inter('reset','inter','linearInterMex');

% inter('reset','inter','splineInterMex');
% configure hyperelasticity based regularization
regularizer('reset','regularizer','mfacid_hysco_hyperEPIstg','alpha',alpha,'beta',beta);

% write parameters to report file
Report('parameters',reportFile,alpha,beta,omega,m,MLdata,minLevel,maxLevel);

% estimate field inhomogeneity based on first image pair
%      (typically acquired without diffusion weighting)
[Bopt,his] = acid_hysco_MLIRepi(MLdata,'minLevel',minLevel,'maxLevel',maxLevel,'plots',plots,'NPIRobj',@acid_hysco_EPINPIRobjFctnStg);

his,
% bring numerical solution to data resolution
%    (may be necessary, if correction is not performed on full data resolution)
Bc  = acid_hysco_prolongate(Bopt,omega,MLdata{maxLevel}.m,m);

% compute transformations and intensity modulations
xc    = reshape(getCellCenteredGrid(omega,m),[],dim);
Bc    = reshape(Bc,m+eye(1,dim));                  % 1-staggered
Bcc   = .5*(Bc(1:end-1,:,:) + Bc(2:end,:,:));   % cell-centered
pB    = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
y1    = xc; y1(:,1) = y1(:,1) + Bcc(:);
y2    = xc; y2(:,1) = y2(:,1) - Bcc(:);
Jac1  = 1 + pB;
Jac2  = 1 - pB;

% save inhomogeneity
fname = acid_bids_filename(V_bu_ref, [keyword '-fmap'], '', '.nii');

VB = V_bu_ref;
VB.descrip = sprintf('HySCOv2: Estimated Imhomogeneity %s',datestr(now()));
VB.dim(phdir) = VB.dim(phdir)+1;
VB.fname = [p_out filesep fname];
VB.dt = [16 0];
VB.private.descrip = VB.descrip;
VB.mat(phdir,4) = VB.mat(phdir,4)-h(phdir)/2;
VB.private.mat = VB.mat;
VB.private.mat0 = VB.mat;
VB.private.dat.dim = VB.dim;
VB.private.dat.fname = VB.fname;

NiB = write_data(reshape(Bc,m+eye(1,dim)), VB, [], [], 1, 1, 'float32', VB.fname);

% apply transformations to original data
I1    = getSplineCoefficients(read_data(V_bu_ref(1)));
I2    = getSplineCoefficients(read_data(V_bd_ref(1)));
I1opt = reshape(splineInterMex(I1,omega,y1).*Jac1(:),m);
I2opt = reshape(splineInterMex(I2,omega,y2).*Jac2(:),m);
Report('images',reportFile,I1,I2,I1opt,I2opt,m)

% get results
rangeB   = [min(Bc(:)), max(Bc(:))];
rangeJac = [min([Jac1(:);Jac2(:)]), max([Jac1(:);Jac2(:)])];
His(1,:) = [100, 100*his.distance(4),his.time,his.iter(minLevel:end),rangeB,rangeJac];
Report(1,reportFile,His(1,:));

% save unwarped blip-up and blip-down image
fname_u2_bu = acid_bids_filename(V_bu_ref, [keyword '-im1'], '_dwi', '.nii');
fname_u2_dw = acid_bids_filename(V_bd_ref, [keyword '-im1'], '_dwi', '.nii');
V_bu_ref.fname = [p_out filesep fname_u2_bu];
V_bd_ref.fname = [p_out filesep fname_u2_dw];
write_data(I1opt, V_bu_ref(1), [], [], 1, 1, 'uint16', V_bu_ref.fname);
write_data(I2opt, V_bd_ref(1), [], [], 1, 1, 'uint16', V_bd_ref.fname);

if ~isempty(V_bu)
    fname_bu = acid_bids_filename(V_bu(1), keyword, '_dwi', '.nii');
    fname_bu = [p_out filesep fname_bu];
end

if ~isempty(V_bd)
    fname_bd = acid_bids_filename(V_bd(1), keyword, '_dwi', '.nii');
    fname_bd = [p_out filesep fname_bd];
end

% now correct remaining image volumes; three cases are covered
% 1) numel(VOI1)==0 and numel(VOI2) > 0      ==> apply correction to VOI2
% 2) numel(VOI1) >0 and numel(VOI2)== 0      ==> apply correction to VOI1
% 3) numel(VOI1)==numel(V0I2) and not(doECC) ==> apply correction to VOI1 and VOI2
% 4) numel(VOI1)==numel(V0I2) and not(doECC) ==> refinement for VOI1 and VOI2

if numel(V_bu)==numel(V_bd)
    for vol = 1:numel(V_bu)
        if doECC
            fprintf('Nonlinear eddy-current correction for volume %d\n',vol);
            
            % read image volumes and normalize intensities
            [I1, I2] = normalizeIntensities(read_data(V_bu(vol)),read_data(V_bd(vol)));
            
            % get multilevel data only for maxLevel
            MLdata = getMultilevel({I1,I2},omega,MLdata{maxLevel}.m,'restrictdim',restrictdim,...
                'minLevel',maxLevel,'maxLevel',maxLevel,'fig',0);
            
            % do iteration only on maxLevel
            [Bdiff,his] = acid_hysco_MLIRepi(MLdata,'minLevel',maxLevel,'maxLevel',maxLevel,'Bc',Bopt,'tolJ',1e-2,...
                'tolG',1e-1,'tolY',1e-1,'plots',plots,'NPIRobj',@acid_hysco_EPINPIRobjFctnStg);
            
            % bring numerical solution to data resolution
            Bdiff = reshape(acid_hysco_prolongate(Bdiff,omega,MLdata{maxLevel}.m,m),m+eye(1,dim));
            
            % update transformation and intensity modulation
            Bdc   = .5*(Bdiff(1:end-1,:,:) + Bdiff(2:end,:,:));   % cell-centered
            pB    = (Bdiff(2:end,:,:) - Bdiff(1:end-1,:,:))/h(1); % partial derivative
            y1    = xc; y1(:,1) = y1(:,1) + Bdc(:);
            y2    = xc; y2(:,1) = y2(:,1) - Bdc(:);
            Jac1  = 1 + pB;
            Jac2  = 1 - pB;

            % save field map
            VBO1 = VB;
            [p_out,fname,~] = fileparts(V_bd(vol).fname);
            VBO1.descrip = sprintf('HySCOv2: Estimated Imhomogeneity %s',datestr(now()));
            VBO1.fname = [p_out filesep fname(1:end-4) '_desc-fmap' '.nii'];
            VBO1.private.descrip = VBO1.descrip;
            VBO1.private.dat.fname = VBO1.fname;
            write_data(reshape(Bdiff,m+eye(1,dim)),VBO1,[],[]);

        else
            fprintf('Apply estimated field inhomogeneity to blip-up and blip-down volumes %d\n',vol);
            Bdiff = Bc;
        end
        
        % load original data
        I1 = getSplineCoefficients(read_data(V_bu(vol)));
        I2 = getSplineCoefficients(read_data(V_bd(vol)));
        
        % apply trafo
        I1opt = reshape(splineInterMex(I1,omega,y1).*Jac1(:),m);
        I2opt = reshape(splineInterMex(I2,omega,y2).*Jac2(:),m);

        % get results
        rangeB   = [min(Bdiff(:)), max(Bdiff(:))];
        rangeJac = [min([Jac1(:);Jac2(:)]), max([Jac1(:);Jac2(:)])];
        if doECC
            His(1+vol,:) = [100*his.distance(3), 100*his.distance(4),his.time,his.iter(minLevel:end),rangeB,rangeJac];
        else
            reduction = SSD(I1opt(:),I2opt(:),omega,m)/SSD(I1(:),I2(:),omega,m);
            His(1+vol,:) = [100, 100*reduction,0,zeros(1,maxLevel-minLevel+1),[-1 1],[-.5 .5]];
        end
        Report(1+vol,reportFile,His(1+vol,:));
        
        if vol==1
            Ni1 = write_data(I1opt, V_bu(vol), VG, [], vol, numel(V_bu), 'uint16', fname_bu);
            Ni2 = write_data(I2opt, V_bd(vol), VG, [], vol, numel(V_bd), 'uint16', fname_bd);           
        else
            Ni1 = write_data(I1opt, V_bu(vol), VG, Ni1, vol, numel(V_bu), 'uint16', fname_bu);
            Ni2 = write_data(I2opt, V_bd(vol), VG, Ni2, vol, numel(V_bu), 'uint16', fname_bd);
        end
    end
    
else
    
    % the number of image volumes does not agree. apply best-known
    % field-inhomogeneity to image volumes
    for vol = 1:numel(V_bu)
        fprintf('Apply estimated field inhomogeneity to blip-up volume %d\n',vol);
        I1    = getSplineCoefficients(read_data(V_bu(vol)));
        I1opt = reshape(splineInterMex(I1,omega,y1).*Jac1(:),m);
        
        if vol==1
            Ni1 = write_data(I1opt, V_bu(vol), VG, [], vol, numel(V_bu), 'uint16', fname_bu);
        else
            Ni1 = write_data(I1opt, V_bu(vol), VG, Ni1, vol, numel(V_bu), 'uint16', fname_bu);
        end
    end

    for vol = 1:numel(V_bd)
        fprintf('Apply estimated field inhomogeneity to blip-down volume %d\n',vol);
        I2 = getSplineCoefficients(read_data(V_bd(vol)));
        I2opt = reshape(splineInterMex(I2, omega,y2).*Jac2(:),m);
        
        if vol == 1
            Ni2 = write_data(I2opt, V_bd(vol), VG, [], vol, numel(V_bd), 'uint16', fname_bd);
        else
            Ni2 = write_data(I2opt, V_bd(vol), VG, Ni2, vol, numel(V_bd), 'uint16', fname_bd);
        end
    end
end

% save json files
if ~isempty(V_bu), acid_save_json(V_bu(1), p_out, keyword); end
if ~isempty(V_bd), acid_save_json(V_bd(1), p_out, keyword); end
    
% save bvals and bvecs files
if ~isempty(V_bu), acid_save_bvals_bvecs(V_bu(1), p_out, keyword); end
if ~isempty(V_bd), acid_save_bvals_bvecs(V_bd(1), p_out, keyword); end

T = toc/60;

T = duration(minutes(T),'format','hh:mm:ss');
disp(['The total time for ' keyword ' was: ' char(T) '.']);


diary off

% finish the report
Report('footer',reportFile);

end

function Dc = SSD(I1, I2, omega, m)

% =========================================================================
% function Dc = SSD(I1,I2,omega,m)
%
% approximates sum-of-squared difference between images I1 and I2 using a
% midpoint rule
%
% Inputs:
%   I1,I2  - image data
%   omega  - computational domain
%   m      - discretization size
%
% Output:
%   Dc     \approx .5 * int_Omega (I1(x)-I2(x))^2 dx
% =========================================================================

    hd = prod((omega(2:2:end)-omega(1:2:end))./m);
    rc = I1-I2;
    Dc = .5*hd*(rc'*rc);

end

function V_tmp = sm_read_vols(V, VG, res)

% =========================================================================
% function sm_write_image(strS,strT,res)
%
% read image volume
%
% Inputs:
%   V    - structure containing image volume information of ith image
%   VG   - structure containing image volume information of target image
%   res  - resampling order
%
% Output:
%   V_tmp - structure containing image volume information of ith image,
%   resliced on the target image
%
% =========================================================================

    V_tmp = acid_read_vols(V, VG, res);
    V_tmp(isnan(V_tmp(:)))=0;
    
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

function [I1,I2] = normalizeIntensities(I1,I2)

% =========================================================================
% function [I1,I2] = normalizeIntensities(I1,I2
%
% normalizes intensities to the range of [0, 256]
%
% Inputs:
%   I1,I2  - image data with intensity range [mini, maxi]
%
% Outputs:
%   I1,I2  - image data with intensity range [0, 256]
% =========================================================================

    mini = min(min(I1(:)),min(I2(:)));
    I1 = I1 - mini;
    I2 = I2 - mini;
    maxi = max(max(I1(:)),max(I2(:)));
    I1 = (256/maxi)*I1;
    I2 = (256/maxi)*I2;

end
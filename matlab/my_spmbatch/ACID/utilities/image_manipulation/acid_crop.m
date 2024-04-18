function [V_out, fname_trafo] = acid_crop(P, matrix, mid_xy, mid_z)

% =========================================================================
% The function cuts the input image to a reduced field of view (rFOV). 
% It will draw a rFOV around the selected point as origin. If edge-length 
% is one- or two-dimensional, it will be cutted to a qubic rFOV. 
% To account for the shift of the origin, the bounding box is calculated
% using a snip of code named world_bb
%
% Inputs:
%   P           - input image
%   matrix      - matrix size of the reduced-FOV
%   mid_xy      - middle points in-plane
%   mid_z       - middle slice
%   p_out       - output folder
%
% Outputs:
%   V_out       - header of the cropped image(s) 
%   fname_trafo - filename of the saved cropping parameters
%
% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
% and Ged Ridgway interpolation tool
% 
% Created by S.Mohammadi 08/10/2014
% Adapted by G.David and B.Fricke
% =========================================================================

%% Check input

% check whether it is the first volume of a 4D series
% if only a single volume was selected
if size(P,1) == 1
    
    tmp = spm_vol([P(1:end-1) '2']);
    
    if ~isempty(tmp)
        struct4D = nifti(P);
        dim4D = struct4D.dat.dim;
        n = dim4D(4);
        if n == 1
            error('A single 3D source image was selected. Choose a single 4D volume or select all 3D volumes manually!');
        end
        P = strcat(repmat(P(1:end-2), n, 1), ',', num2str([1:n]'));
    end
end

% load in header(s) and image(s)
V  = spm_vol(P);
I  = spm_read_vols(V);

% flip if necessary
Sign = sign([V(1).mat(1,1),V(1).mat(2,2),V(1).mat(3,3)]);
if Sign(1) == -1, I = flip(I,1); end
if Sign(2) == -1, I = flip(I,2); end
if Sign(3) == -1, I = flip(I,3); end

% get first image
I1 = I(:,:,:,1);

% check feasibility of cropping
if (matrix(1)<3 && matrix(1)~=-1) || (matrix(2)<3 && matrix(2)~=-1)
    error('Incorrect matrix dimensions: the matrix size in the x and y dimensions has to be at least 3.');
elseif (matrix(3)<1 && matrix(3)~=-1)
    error('Incorrect matrix dimensions: the matrix size in the z dimension has to be at least 1.');    
end

% define output directory
keyword = 'CROP';
[path,fname,~] = spm_fileparts(V(1).fname);

p_out = acid_bids(path,fname,keyword,1);


%% Determine center of the image

% z coordinate
if ~numel(mid_z)==0
    if (mid_z(1) < 1 || mid_z(1) > size(I1,3))
        error('Incorrect midpoint specification: the entered slice does not lie in the image.');
    end
else
     % if not defined, select the middle slice
    mid_z = round(size(I1,3)/2);
end

% x and y coordinates
if ~numel(mid_xy)==0
    mid_xy = round(mid_xy);
    if (mid_xy(1)<1 || mid_xy(1)>size(I1,1))
        error('Incorrect midpoint specification: the midpoint in x does not lie in the image.');
    elseif (mid_xy(2)<1 || mid_xy(2)>size(I1,2))
        error('Incorrect midpoint specification: the midpoint in y does not lie in the image.');
    end  
else
    % if not defined, select in figure
    Fig_1 = figure('Name','Cropping'); colormap gray;
    axis square, title('Select center position!');
    imshow(rot90(squeeze(I1(:,:,mid_z))),[0 max(I1(:))]);
    [cent_x, cent_y] = ginput(1);
    mid_xy(1) = round(cent_x);
    mid_xy(2) = round(size(I1,2)-cent_y+1); % flip because indexing starts from the bottom of the image, while ginput starts from the top
    close(Fig_1)
end

%% Determine the FOV of the cropped images
% x-FOV
if matrix(1) == -1
    xmin = 1; xmax = size(I1,1);
else
    xmin = mid_xy(1) - ceil(matrix(1)/2) + 1; 
    if xmin < 1, xmin = 1; end
    if xmin > size(I1,1), xmin = size(I1,1); end
    xmax = mid_xy(1) + floor(matrix(1)/2);
    if xmax < 1, xmax = 1; end
    if xmax > size(I1,1), xmax = size(I1,1); end
end

% y-FOV
if matrix(2) == -1
    ymin = 1; ymax = size(I1,2);
else
    ymin = mid_xy(2) - ceil(matrix(2)/2) + 1;
    if ymin < 1, ymin = 1; end
    if ymin > size(I1,2), ymin = size(I1,2); end
    ymax = mid_xy(2) + floor(matrix(2)/2);
    if ymax < 1, ymax = 1; end
    if ymax > size(I1,2), ymax = size(I1,2); end
end

% z-FOV
if matrix(3) == -1
    zmin = 1; zmax = size(I1,3);
else
    zmin = mid_z - ceil(matrix(3)/2) + 1;
    if zmin < 1, zmin = 1; end
    if zmin > size(I1,3), zmin = size(I1,3); end
    zmax = mid_z + floor(matrix(3)/2);
    if zmax < 1, zmax = 1; end
    if zmax > size(I1,3), zmax = size(I1,3); end
end

% save matrix size and midpoint
fname = acid_bids_filename(V(1), keyword);
fname_trafo = [p_out filesep fname(1:end-4) '-params' '.mat'];
save(fname_trafo,'matrix','mid_xy','mid_z');

%% Create and write out cropped images
if Sign(1) == -1, a = size(I1,1)-xmax+1; else, a = xmin; end
if Sign(2) == -1, b = size(I1,2)-ymax+1; else, b = ymin; end
if Sign(3) == -1, c = size(I1,3)-zmax+1; else, c = zmin; end
 
% initialization
I_crop = zeros(length(xmin:xmax),length(ymin:ymax),length(zmin:zmax),size(P,1));

% modify header for the cropped image(s)
VG = V(1);
[p,f,e]  = fileparts(V(1).fname);
VG.fname = fullfile(p,[f e]);
VG.dim(1:3) = [numel(xmin:xmax) numel(ymin:ymax) numel(zmin:zmax)];
VG.mat      = V(1).mat*spm_matrix([a-1,b-1,c-1]);

% crop image(s)
for i = 1:size(P,1)
    I_crop(:,:,:,i) = I(xmin:xmax,ymin:ymax,zmin:zmax,i);    
end

% flip dimension(s) if necessary
if Sign(1) == -1, I_crop = flip(I_crop,1); end
if Sign(2) == -1, I_crop = flip(I_crop,2); end
if Sign(3) == -1, I_crop = flip(I_crop,3); end

%  cropped image
V_out = acid_write_vol(I_crop, VG, p_out, 'CROP', 'same','','',1);

% save json file
acid_save_json(VG, p_out, keyword);

% save bvals and bvecs
acid_save_bvals_bvecs(VG, p_out, keyword);

disp('Done!');

end
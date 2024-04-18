function Vout = acid_write_vol(varargin)

% =========================================================================
% The script writes out a 3D or 4D nifti file received in the input (I) using the header information of V.
%
% Usage:  Vout = acid_write_vol(I, V, p_out, keyword, suffix, mask, desc, dummy_float)
%
% Inputs:
%   I           - 3D or 4D matrix, to be written as a nifti file
%   V           - reference header used to create the header of the output nifti file
%   p_out       - folder for output files
%   keyword     - tag for the 'desc' field of the BIDS output filename (e.g., sub-*_desc-ECMOCO_*.nii)
%   suffix      - postfix for the BIDS outout filename (appended to the input filename)
%   mask        - mask applied on I
%   desc        - description entry of the header
%   dummy_float - dummy for floating number precision (0 for uint16, 1 for float32)
%
% Output:
%   Vout    - header of the output nifti file
%   Plus, a nifti file is saved with the specified name.
%
% created by S. Mohammadi
% adapted by G. David, B. Fricke
% =========================================================================

% get number of supplied arguments
[~,nargin] = size(varargin);

% input #1 (obligatory)
I = varargin{1};

% input #2 (obligatory)
V = varargin{2};

% create nifti object 
Ni = nifti;
Ni.mat  = V(1).mat;
Ni.mat0 = V(1).mat;

% input #3 (obligatory)
if nargin >= 3 && ~isempty(varargin{3})
    p_out = varargin{3};
end

% input #4
if nargin >= 4 && ~isempty(varargin{4})
    keyword = varargin{4};      
else
    keyword = '';   
end

% input #5
if nargin >= 5 && ~isempty(varargin{5})
    suffix = varargin{5};
    if strcmp(suffix, 'same')
        pos1 = find(V.fname == '_', 1, 'last');
        pos2 = find(V.fname == '.', 1, 'last');
        tmp = V.fname(pos1:pos2-1);
        if max(contains(['_dwi','_map','_T1','_T2'], tmp))
            suffix = tmp;
        else
            suffix = '';
        end
    end
else
    suffix = '';   
end

% input #6
if nargin >= 6 && ~isempty(varargin{6})
    mask       = varargin{6};
    vol1       = zeros(V.dim);
    vol1(mask) = I;       
    I          = vol1;
end

% input #7
if nargin >= 7 && ~isempty(varargin{7})
    Ni.descrip = varargin{7};
end

% input #8
if nargin >= 8 && ~isempty(varargin{8})
    dummy_float = varargin{8};
else
    dummy_float = 0;
end
if dummy_float == 0
    dt = [spm_type('uint16'),spm_platform('bigend')];
else
    dt = [spm_type('float32'),spm_platform('bigend')];
end

% generate filename
fname = acid_bids_filename(V, keyword, suffix, '.nii');
fname = [p_out filesep fname];

% create output nifti file
Ni.dat = file_array(fname,size(I),dt,0,1,0);
create(Ni);
if numel(size(I))==3
    spm_progress_bar('Init',size(I,3),Ni.descrip,'planes completed');
    for p = 1:size(I,3)
        Ni.dat(:,:,p) = I(:,:,p);
        spm_progress_bar('Set',p);
    end
elseif numel(size(I))==4
    spm_progress_bar('Init',size(I,4),Ni.descrip,'volumes completed');
    for p = 1:size(I,4)
        Ni.dat(:,:,:,p) = I(:,:,:,p);
        spm_progress_bar('Set',p);
    end
end

spm_progress_bar('Clear');
Vout = spm_vol(fname);

end
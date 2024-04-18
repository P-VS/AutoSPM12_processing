function fname = acid_bids_filename(varargin)

% Usage: acid_bids_filename(V, keyword, suffix, ext)

% input #1 (obligatory)
V = varargin{1};
V = V(1);

% input #2 (obligatory)
if nargin >= 2 && ~isempty(varargin{2})
    keyword = varargin{2};      
else
    keyword = '';   
end

% input #3
if nargin >= 3 && ~isempty(varargin{3})
    suffix = varargin{3};
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

% input #4
if nargin >= 4 && ~isempty(varargin{4})
    ext = varargin{4};      
else
    ext = '.nii';   
end

% remove file extension
[~, fname, ~] = spm_fileparts(V(1).fname);

% determine whether there is suffix
if contains(fname, '_dwi')
    dummy_suffix  = 1;
    length_suffix = 4;
elseif contains(fname, '_map')
    dummy_suffix  = 1;
    length_suffix = 4;
elseif contains(fname, '_T1')
    dummy_suffix  = 1;
    length_suffix = 3;
elseif contains(fname, '_T2')
    dummy_suffix  = 1;
    length_suffix = 3;
else
    dummy_suffix = 0;
end

% desc field, suffix
if contains(fname, '_desc-') && dummy_suffix && size(fname,2) > 8
    if ~isempty(keyword)
        fname = [fname(1:end-length_suffix) '-' keyword suffix ext];
    else
        fname = [fname(1:end-length_suffix) keyword suffix ext];
    end  
% desc field, no suffix 
elseif contains(fname, '_desc-') && ~dummy_suffix
    if ~isempty(keyword)
        fname = [fname '-' keyword suffix ext];
    else
        fname = [fname suffix ext];
    end
% no desc field, suffix
elseif ~contains(fname, '_desc-') && dummy_suffix
    fname = [fname(1:end-length_suffix) '_desc-' keyword suffix ext];  
% no desc field, no suffix
else
    fname = [fname '_desc-' keyword suffix ext];
end
end
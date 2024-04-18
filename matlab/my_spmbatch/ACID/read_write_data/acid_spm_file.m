function varout = acid_spm_file(varargin)

if(nargin>=1)
    filename=char(varargin{1});
    for i=1:size(filename,1)
        if(nargin>=3)
            options=varargin{2};
            prename=varargin{3};
            if(strcmp('prefix',options))
                [p,n,e] = spm_fileparts(filename(i,:));
                varout{i}  = fullfile(p,[prename, n, e]);
                if(nargin>=4 && nargin<=5)
                    options=varargin{4};
                    endname=varargin{5};
                    if(strcmp('format',options))
                        [p,n,e]=spm_fileparts(filename(i,:));
                        varout{i}  = fullfile(p,[prename,n, endname]);
                    elseif(strcmp('ending',options))
                        varout{i}  = fullfile(p,[prename,n, endname, e]);
                    end
                elseif(nargin>=7)
                    options=varargin{4};
                    endname=varargin{5};
                    options2=varargin{6};
                    endname2=varargin{7};
                    if(strcmp('ending',options))
                        if(strcmp('format',options2))
                            varout{i}  = fullfile(p,[prename,n, endname, endname2]);
                        else
                            error('Error in assigning dependencies');
                        end
                    else
                        error('Error in assigning dependencies');
                    end
                end
            elseif(strcmp('ending',options))
                [p,n,e]=spm_fileparts(filename(i,:));
                varout{i}  = fullfile(p,[n prename e]);
                if(nargin>=5)
                    options=varargin{4};
                    endname=varargin{5};
                    if(strcmp('format',options))
                        [p,n,e]=spm_fileparts(filename(i,:));
                        varout{i}  = fullfile(p,[n prename endname]);
                    end
                end
            end
        else
            varout  = varargin{1};
        end
    end
end
end
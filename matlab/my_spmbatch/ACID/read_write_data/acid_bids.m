function p_out = acid_bids(path, fname, keyword, dummy_derivatives)
% This function generates a BIDS compliant folder structure for the ACID Toolbox.
% B.Fricke 11.07.2022

cd(path);

dummy_bidsname = 0;
dummy_bids_storage = 0;

if contains(fname,"sub-")
    dummy_bidsname = 1;
    fname_split = split(fname,"_");
    
    if strcmp(fname_split{1}(1:4),"sub-")
        subject_folder_name = fname_split{1};
    end
end

if contains(path,"sub-")
    if dummy_bidsname == 1
        dummy_bids_storage = 1;
    end

    path_split = split(path,filesep);
    
    i = 0;
    
    while contains(path_split{end-i},"sub-") == 0
        i = i + 1;
    end
    i = i + 1;
end


cd('../')
path_possible_derivatives_folder = cd;
% cd('../')
% cd('../')
% path_possible_derivatives_bids_dataset_folder = cd;
cd(path);
if strcmp(path_possible_derivatives_folder(end-10:end), 'derivatives')
    cd('../');
    path_derivatives = cd;
    path = path_derivatives;
    dummy = 1;

% elseif strcmp(path_possible_derivatives_bids_dataset_folder(end-10:end), 'derivatives')
%     cd('../');
%     path_derivatives = cd;
%     path = path_derivatives;
%     dummy = 1;

elseif dummy_bids_storage
    for counter = 1:i
        cd('../');
    end
    path_possible_derivatives_bids_dataset_folder = cd;
    path = [path_possible_derivatives_bids_dataset_folder filesep 'derivatives'];
    if not(isfolder(path)), mkdir(path); end
    path = [path_possible_derivatives_bids_dataset_folder filesep 'derivatives' filesep 'ACID'];
    if not(isfolder(path)), mkdir(path); end   
    path = [path_possible_derivatives_bids_dataset_folder filesep 'derivatives' filesep 'ACID' filesep subject_folder_name];
    if not(isfolder(path)), mkdir(path); end   
    dummy = 2;
    
elseif dummy_derivatives
    path_bids_higher_higher_level = cd;
    path = [path_bids_higher_higher_level filesep 'derivatives'];
    if not(isfolder(path)), mkdir(path); end
    dummy = 2;
    
else
    p_out = path;
    dummy = 3;
    
end

% 
%   messung/dwi/sub1-ecomo.nii
% 


if dummy==1 || dummy==2
    p_out = [path filesep keyword];
    if not(isfolder(p_out))
        mkdir(p_out)
    else
        counter = 2;
        p_out = [path filesep keyword '-Run_' num2str(counter)];

        while isfolder(p_out)
            counter = counter + 1;
            p_out = [path filesep keyword '-Run_' num2str(counter)];
        end
        mkdir(p_out)
    end
end
end
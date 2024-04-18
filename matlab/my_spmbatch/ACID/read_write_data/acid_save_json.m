function acid_save_json(V, p_out, keyword)

    % get filename of json associated with the input files
    [p, fname, ~] = spm_fileparts(V(1).fname);
    fname_json = [p filesep fname '.json'];
    
    % generate new filenames for bvals and bvecs
    fname_json_new = acid_bids_filename(V, keyword, '_dwi', '.json');
    fname_json_new = [p_out filesep fname_json_new];

    % save json file
    try decoded_json = fileread(fname_json);

        encoded_json = jsondecode(decoded_json);
        if ~isfield(encoded_json,'Description')
            encoded_json.Description = [];
        end

        if ~isfield(encoded_json,'Sources')
            encoded_json.Description = [];
        end

        encoded_json.Sources{1} = V(1).fname;
        encoded_json.Description{end+1} = keyword;
        spm_jsonwrite(fname_json_new, encoded_json, 'indent', '\t');
    catch
    end
end
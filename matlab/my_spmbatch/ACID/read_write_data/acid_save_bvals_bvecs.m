function acid_save_bvals_bvecs(V, p_out, keyword)

    % get filenames of bvals and bvecs associated with the input files
    [p, fname, ~] = spm_fileparts(V(1).fname);
    fname_bvals = [p filesep fname '.bval'];
    fname_bvecs = [p filesep fname '.bvec'];

    % generate new filenames for bvals and bvecs
    fname_bvals_new = acid_bids_filename(V, keyword, '_dwi', '.bval');
    fname_bvecs_new = acid_bids_filename(V, keyword, '_dwi', '.bvec');

    % save bvals and bvecs under the new names
    try bvals = load(fname_bvals);
        if size(bvals,1)>1
            bvals = bvals';
        end
        dlmwrite([p_out filesep fname_bvals_new], bvals);
    catch
    end

    try bvecs = load(fname_bvecs);
        dlmwrite([p_out filesep fname_bvecs_new], bvecs);
    catch
    end

end
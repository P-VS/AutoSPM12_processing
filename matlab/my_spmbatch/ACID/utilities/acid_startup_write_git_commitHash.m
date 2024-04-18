function acid_startup_write_git_commitHash(p_out)

    path_git_repo = fileparts(mfilename('fullpath'));
    % path = [p_out filesep];
    
    cd(path_git_repo);
    
    [status, commitHash] = system('git rev-parse HEAD');

    cd(p_out)

    if status ~= 0
        disp('Unable to get commit hash');
    else
        disp(['Current commit hash is: ' commitHash]);
    end

end

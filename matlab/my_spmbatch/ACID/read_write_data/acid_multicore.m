function npool = acid_multicore(npool)
% B.Fricke 06.02.2024


if npool>1
    try p = gcp('nocreate');
    if p.NumWorkers ~= npool
        delete(gcp('nocreate'));
        if exist('parpool')>0
            poolobj = parpool('local',npool);
            dummy_matlabpool = false;
        elseif exist('matlab')>0
            parpool('open',npool);
            dummy_matlabpool = true;
        else
            npool = 1;
            warning('No parallel processing possible on this matlab version! Proceed without parallel processing.')
        end
    end
    catch
        % delete(gcp('nocreate'));
        if exist('parpool')>0
            poolobj = parpool('local',npool);
            dummy_matlabpool = false;
        elseif exist('matlab')>0
            parpool('open',npool);
            dummy_matlabpool = true;
        else
            npool = 1;
            warning('No parallel processing possible on this matlab version! Proceed without parallel processing.')
        end
    end

elseif npool == 0
    try p = gcp('nocreate');
        npool = p.NumWorkers;
    catch
        warning('No parallel pool was found! Proceed without parallel processing.')
    end

end


end
function my_spmbatch_noparallel(sublist,nsessions,task,datpath,params,save_intermediate_results)

if params.do_bpfilter || params.do_aCompCor || params.do_ICA_AROMA || params.do_noiseregression || params.do_DENN
    do_denoising = true;

    if ~params.do_onlydenoise
        params.prefix = 'e';
        if params.fieldmap, params.prefix = ['u' params.prefix]; end
        if params.do_realignment && ~params.fieldmap, params.prefix = ['r' params.prefix]; end
        if params.pepolar, params.prefix = ['u' params.prefix]; end
        if params.do_slicetime, params.prefix = ['a' params.prefix]; end
        if params.meepi && ~contains(params.combination,'none')
            params.prefix = ['c' params.prefix]; 
            params.mecombined = true;
        else
            params.mecombined = false;
        end
        if params.do_normalization, params.prefix = ['w' params.prefix]; end
    end
else
    do_denoising = false;
end

for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        for k=1:numel(task)
            itstart = tic;

            %% make and run batch

            if ~params.do_onlydenoise
                [delfiles,keepfiles] = my_spmbatch(sublist(i),nsessions(j),task{k},datpath,params);
            
                %% Clean up unnecessary files
                cleanup_intermediate_files(sublist(i),nsessions(j),datpath,delfiles,keepfiles,save_intermediate_results,params);

                %% Print and save realignment paramers  
                save_rp_plot(sublist(i),nsessions(j),task{k},datpath,params);
            end

            if do_denoising
                [delfiles,keepfiles] = my_spmbatch_onlydenoise(sublist(i),nsessions(j),task{k},datpath,params);

                %% Clean up unnecessary files
                cleanup_intermediate_files(sublist(i),nsessions(j),datpath,delfiles,keepfiles,save_intermediate_results,params);
            end

            itstop = toc(itstart);

            fprintf(['subject ' num2str(sublist(i)) ' session ' num2str(nsessions(j)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
        end
    end
end
function my_spmbatch_parallel(sublist,nsessions,task,datpath,params,save_intermediate_results)

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

datlist = zeros(numel(sublist)*numel(nsessions),2);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        datlist(dpos,1) = sublist(i);
        datlist(dpos,2) = nsessions(j);

        dpos = dpos+1;
    end
end

for k = 1:numel(task)
    pa=parpool(min([3,numel(datlist(:,1))])); %25 is the maximum number of workers allowed in the 'local' profile while 10 is set to avoid memory issues on my computer
    parfor i = 1:numel(datlist(:,1))
        %try
            if ~params.do_onlydenoise
                %% make and run batch
                [delfiles,keepfiles] = my_spmbatch(datlist(i,1),datlist(i,2),task{k},datpath,params);
                
                %% Clean up unnecessary files
                cleanup_intermediate_files(datlist(i,1),datlist(i,2),datpath,delfiles,keepfiles,save_intermediate_results,params);  
            end

            %% Do denoising
            if do_denoising

                [delfiles,keepfiles] = my_spmbatch_onlydenoise(datlist(i,1),datlist(i,2),task{k},datpath,params);

                %% Clean up unnecessary files
                cleanup_intermediate_files(datlist(i,1),datlist(i,2),datpath,delfiles,keepfiles,save_intermediate_results,params);
            end
        %catch
         %   fprintf(['\nError for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' task ' task{k} '\n'])
        %end
    end
    delete(pa)

    if ~params.do_onlydenoise
        for i = 1:numel(datlist(:,1))
            %% Print and save realignment paramers  
            save_rp_plot(datlist(i,1),datlist(i,2),task{k},datpath,params);
        end
    end
end
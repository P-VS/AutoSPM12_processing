function my_spmbatch_parallel(sublist,nsessions,task,datpath,params,save_intermediate_results)

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
        try
            %% make and run batch
            
            if ~params.do_onlydenoise
                [delfiles,keepfiles] = my_spmbatch(datlist(i,1),datlist(i,2),task{k},datpath,params);
                
                %% Clean up unnecessary files
                cleanup_intermediate_files(datlist(i,1),datlist(i,2),datpath,delfiles,keepfiles,save_intermediate_results,params);
        
            else
                my_spmbatch_onlydenoise(sublist(i),nsessions(j),task{k},datpath,params);
            end
        catch
            fprintf(['\nError for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' task ' task{k} '\n'])
        end
    end
    delete(pa)

    for i = 1:numel(datlist(:,1))
        %% Print and save realignment paramers  
        save_rp_plot(datlist(i,1),datlist(i,2),task{k},datpath,params);
    end
end
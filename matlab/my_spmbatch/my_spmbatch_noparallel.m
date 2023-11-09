function my_spmbatch_noparallel(sublist,nsessions,task,datpath,params,save_intermediate_results)

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
            else
                my_spmbatch_onlydenoise(sublist(i),nsessions(j),task{k},datpath,params);
            end

            itstop = toc(itstart);

            fprintf(['subject ' num2str(sublist(i)) ' session ' num2str(nsessions(j)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
        end
    end
end
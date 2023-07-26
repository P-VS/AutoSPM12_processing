function my_spmbatch_parallelprocessing(sublist,nsessions,task,params)

datlist = zeros(numel(sublist)*numel(nsessions),2);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        datlist(dpos,1) = sublist(i);
        datlist(dpos,2) = nsessions(j);

        dpos = dpos+1;
    end
end

pa=parpool(min([3,numel(datlist(:,1))])); %25 is the maximum number of workers allowed in the 'local' profile while 10 is set to avoid memory issues on my computer
parfor i = 1:numel(datlist(:,1))
            
    %% make batch
    matlabbatch = my_spmprocessingbatch(datlist(i,1),datlist(i,2),task,params);

    spm_jobman('run', matlabbatch);
end
delete(pa)
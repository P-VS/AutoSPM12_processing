function my_spmbatch_noparallelprocessing(sublist,nsessions,task,params)

for i = sublist
    for j = nsessions
        matlabbatch = my_spmprocessingbatch(i,j,task,params);

        spm_jobman('run', matlabbatch);
    end
end
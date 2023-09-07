function my_spmdcm2niibatch(sub,params,useparfor)

%Code based on dcm2nii from
%(https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)
%(https://www.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii)

substring = ['sub-' num2str(sub,'%02d')];

for si=1:numel(params.mridata)
    subpath = fullfile(params.datpath,substring,['ses-' num2str(params.mridata(si).session,'%03d')]);
    
    infolder = fullfile(params.datpath,substring,params.mridata(si).folder);
    outfolder = fullfile(subpath,params.mridata(si).seqtype);

    if ~isfolder(outfolder)
        mkdir(outfolder);
    end

    %% nii file name
    outfname = [substring '_' params.mridata(si).name];

    my_spmdicm2niix(infolder, outfolder, '.nii', outfname, useparfor)
    
end
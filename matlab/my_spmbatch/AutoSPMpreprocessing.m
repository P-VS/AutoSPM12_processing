function AutoSPMpreprocessing

%Script to do the auto preprocessing in SPM12
%
%Preparation:
% Convert the DICOM files into nifti using dcm2niix in MROCroGL
% For the anatomical scans, set 'Crop 3D Images' on
%
%* Organise the data in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1)
%                   Files: sub-##_T1w.nii and sub-##_T1w.json
%                -func: containes the fmri data
%                   Files: sub-##_task-..._bold.nii and sub-##_task-..._bold.json
%                -fmap: containnes the gradient pololarity (blip-up/down) filpt data or the fieldmap scans
%                   Files in case of inverted gradient polarity: sub-##_dir-pi_epi.nii and sub-##_dir-pi_epi.json
%                   Files in case of fieldmap scans: (image 1 in file is amplitude, image 2 in file is phase)
%                          sub-##_fmap_echo-1.nii and sub-##_fmap_echo-1.json
%                          sub-##_fmap_echo-2.nii and sub-##_fmap_echo-2.json
    
%* IMPORTANT: !! Look at your data before starting any (pre)processing. Losing time in trying to process bad data makes no sense !!

%Script written by dr. Peter Van Schuerbeek (Radiology UZ Brussel)

warnstate = warning;
warning off;

%spm_figure('Clear','Interactive')

% User interface.
SPMid                 = spm('FnBanner',mfilename,'2.10');
[Finter,Graf,CmdLine] = spm('FnUIsetup','Preproces SPM');

%%
%%Give the basic input information of your data

datpath = '/Volumes/LaCie/UZ_Brussel/Labo_fMRI/Full_dataset/';

first_sub = 7;
last_sub = 7;
sublist = [2];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0

task ={'affect_run-1'};

params.nechoes = 1; %number of echoes for ME-fMRI. The combination is done wi=TEi 
params.dummytime = 10; %time in seconds

params.reorient = true;

%fieldmap or pepolar should be true, the other should be false
params.fieldmap = false;
params.pepolar = false;

params.do_realignment = true;
params.do_slicetime = true;
params.do_segmentation = false;
params.do_normalization = true;
params.do_mot_derivatives = true;
params.do_aCompCor = true;
params.do_noiseregression = true;
params.do_ICA_AROMA = true;
params.do_smoothing = true;

params.normvox = [1.5 1.5 1.5];
params.smoothfwhm = 6;

use_parallel = false;
save_intermediate_results = false;

%% BE CAREFUL WITH CHANGING THE CODE BELOW THIS LINE !!
%---------------------------------------------------------------------------------------
fprintf('Start with preprocessing \n')

curdir = pwd;

spm('defaults', 'FMRI');

if use_parallel
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
        pa=parpool(min([10,numel(datlist(:,1))])); %25 is the maximum number of workers allowed in the 'local' profile while 10 is set to avoid memory issues on my computer
        parfor i = 1:numel(datlist(:,1))
            try
                %% make and run batch
                
                [delfiles,keepfiles] = my_spmbatch(datlist(i,1),datlist(i,2),task{k},datpath,params);
                
                %% Clean up unnecessary files
                cleanup_intermediate_files(datlist(i,1),datlist(i,2),datpath,delfiles,keepfiles,save_intermediate_results);
        
            catch
                fprintf(['\nError for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' task ' task{k} '\n'])
            end
        end
        delete(pa)

        for i = 1:numel(datlist(:,1))
            %% Print and save realignment paramers  
            save_rp_plot(datlist(i,1),datlist(i,2),task{k},datpath);
        end
    end
else
    for i = 1:numel(sublist)
        for j = 1:numel(nsessions)
            for k=1:numel(task)
                itstart = tic;

                %% make and run batch

                [delfiles,keepfiles] = my_spmbatch(sublist(i),nsessions(j),task{k},datpath,params,save_intermediate_results);

                %% Clean up unnecessary files
                cleanup_intermediate_files(sublist(i),nsessions(j),datpath,delfiles,keepfiles,save_intermediate_results);

                %% Print and save realignment paramers  
                save_rp_plot(sublist(i),nsessions(j),task{k},datpath);

                itstop = toc(itstart);

                fprintf(['subject ' num2str(sublist(i)) ' session ' num2str(nsessions(j)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
            end
        end
    end
end

fprintf(['\nDone\n'])

end
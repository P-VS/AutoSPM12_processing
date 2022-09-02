function my_spmbatch_normfslsegments

%Script to do the auto normalizze the segmentation maps found with FSL
%
%Preparation:
% Convert the DICOM files into nifti using dcm2niix in MROCroGL
% For the anatomical scans, set 'Crop 3D Images' on
%
%* Organise the data in BIDS format
%    - datpath
%        -sub-##
%            -ses-00# (if your experiment contains multiple session per subject)
%                -anat: containes the anatomical data (3D T1) and
%                segmentation maps
%                   Files: sub-##_T1w.nii and sub-##_T1w.json
%
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

datpath = '/Volumes/LaCie/UZ_Brussel/DeNN_motor/Data';

first_sub = 7;
last_sub = 7;
sublist = [2:5];%list with subject id of those to preprocess separated by , (e.g. [1,2,3,4]) or alternatively use sublist = [first_sub:1:last_sub]
nsessions = [1]; %nsessions>0

params.reorient = true;

params.normvox = [1.5 1.5 1.5];

use_parallel = true;
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

    pa=parpool(min([10,numel(datlist(:,1))])); %25 is the maximum number of workers allowed in the 'local' profile while 10 is set to avoid memory issues on my computer
    parfor i = 1:numel(datlist(:,1))
        try
            %% make and run batch
            [delfiles,keepfiles] = my_spmbatch_normfslsegments_job(datlist(i,1),datlist(i,2),datpath,params);

            %% Clean up unnecessary files
            cleanup_intermediate_files(datlist(i,1),datlist(i,2),datpath,delfiles,keepfiles,save_intermediate_results);
        catch
            fprintf(['\nError for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) '\n'])
        end
    end
    delete(pa)
else
    for i = 1:numel(sublist)
        for j = 1:numel(nsessions)
            itstart = tic;

            %% make and run batch
            [delfiles,keepfiles] = my_spmbatch_normfslsegments_job(sublist(i),nsessions(j),datpath,params);

            %% Clean up unnecessary files
            cleanup_intermediate_files(sublist(i),nsessions(j),datpath,delfiles,keepfiles,save_intermediate_results);

            itstop = toc(itstart);

            fprintf(['subject ' num2str(sublist(i)) ' session ' num2str(nsessions(j)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
        end
    end
end

fprintf(['\nDone\n'])

end
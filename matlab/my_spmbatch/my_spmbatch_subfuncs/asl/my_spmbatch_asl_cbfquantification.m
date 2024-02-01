function [ppparams,delfiles,keepfiles] = my_spmbatch_asl_cbfquantification(ppparams,delfiles,keepfiles)

%% Step 1: make head mask

M0I = spm_vol(ppparams.subm0scan);
m0dat = spm_read_vols(M0I);

m0mask = my_spmbatch_mask(m0dat);

% from Alsop 2015 MRM
BloodT1 = 1.650; 
LabelingEfficiency = 0.85; %labeling efficiency for PCASL
SupressionEfficiency = 0.75; %Effect of background suppression on labeled spins (0.75 for GE 3D PCASL)
Lambda = 0.9; %blood-brain partition coefficient
SCF = 32; %@GE the deltam images are upscalled by a factor 32

jsondat = fileread(ppparams.subdmjson);
jsondat = jsondecode(jsondat);

LabelingDuration = jsondat.LabelingDuration;
PLD = jsondat.PostLabelDelay;
NumberOfAverages = jsondat.NumberOfAverages;

%% Step 2: Calculating the CBF images

% CBF = (6000 * Lambda * DM * exp(PLD/BloodT1)) ./ (2 * LabelingEfficiency * SupressionEfficiency * SCF * NumberOfAverages * BloofT1 * M0 * (1-exp(-LabelingDuration/BloodT1)))

DM = spm_vol(ppparams.subdeltam);
dmdat = double(spm_read_vols(DM));

M0 = spm_vol(ppparams.subm0scan);
m0dat = double(spm_read_vols(M0));
%m0dat = m0dat * M0.pinfo(1,1) + M0.pinfo(2,1);

cbfdat = (6000*Lambda*dmdat*exp(PLD/BloodT1)) ./ (2*LabelingEfficiency*SupressionEfficiency*SCF*NumberOfAverages*BloodT1*m0dat*(1-exp(-LabelingDuration/BloodT1)));
cbfdat = cbfdat .* (m0mask>0);

[dmpth,dmname,~] = fileparts(ppparams.subdeltam);
dmparts = split(dmname,'deltam');

cbffile = fullfile(dmpth,[dmparts{1} 'cbf.nii']);

CBF = DM;
CBF.fname = cbffile;
CBF.descrip = 'CBF';
CBF = rmfield(CBF,'pinfo');
CBF = spm_create_vol(CBF);
CBF = spm_write_vol(CBF,cbfdat);

ppparams.subcbf = cbffile;
delfiles{numel(delfiles)+1} = {cbffile};
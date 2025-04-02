function [ppparams,delfiles,keepfiles] = my_spmbatch_split_asl_bold(params,ppparams,ie,delfiles,keepfiles)

fprintf('Split ASL/BOLD \n')

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
funcdat = spm_read_vols(Vfunc);

jsondat = fileread(ppparams.func(ppparams.echoes(1)).jsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;
Ny = 1/(2*tr);

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

[bolddat,~] = fmri_cleaning(funcdat(:,:),0,[tr 0.008 Ny-0.008],[],[],'restoremean','on');
bolddat = reshape(bolddat(:,:),s);

fname = split(ppparams.func(ie).funcfile,'_aslbold.nii');

Vbold = Vfunc;
for iv=1:numel(Vbold)
    Vbold(iv).fname = fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_aslbold.nii']);
    Vbold(iv).descrip = 'my_spmbatch - split bold';
    Vbold(iv).pinfo = [1,0,0];
    Vbold(iv).n = [iv 1];
end

Vbold = myspm_write_vol_4d(Vbold,bolddat);

[nbolddat,~] = fmri_cleaning(funcdat(:,:),0,[tr Ny-0.008 Ny],[],[],'restoremean','on');
nbolddat = reshape(nbolddat(:,:),s);

Vasl = Vfunc;
for iv=1:numel(Vasl)
    Vasl(iv).fname = fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_label.nii']);
    Vasl(iv).descrip = 'my_spmbatch - split asl';
    Vasl(iv).pinfo = [1,0,0];
    Vasl(iv).n = [iv 1];
end

Vasl = myspm_write_vol_4d(Vasl,nbolddat);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_aslbold.nii'])};
delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_label.nii'])};

ppparams.func(ie).funcfile = [fname{1} '_aslbold.nii'];
ppparams.func(ie).prefix = ['f' ppparams.func(ie).prefix];

ppparams.func(ie).perffile = [ppparams.func(ie).prefix fname{1} '_label.nii'];

clear funcdat bolddat Vfunc Vbold Vasl
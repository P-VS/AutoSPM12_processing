function [ppparams,delfiles,keepfiles] = my_spmbatch_split_asl_bold(params,ppparams,ie,delfiles,keepfiles)

fprintf('Split ASL/BOLD \n')

if ~exist(fullfile(ppparams.subpath,'perf'),'dir'), mkdir(fullfile(ppparams.subpath,'perf')); end

ppparams.subperfdir = fullfile(ppparams.subpath,'perf');

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
funcdat = spm_read_vols(Vfunc);

jsondat = fileread(ppparams.func(ppparams.echoes(1)).jsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

[bolddat,~] = fmri_cleaning(funcdat(:,:),0,[tr 0 0.1],[],[],'restoremean','on');

bolddat = reshape(bolddat(:,:),s);

fname = split(ppparams.func(ie).funcfile,'_aslbold.nii');

Vbold = Vfunc;
for iv=1:numel(Vbold)
    Vbold(iv).fname = fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_bold.nii']);
    Vbold(iv).descrip = 'my_spmbatch - split bold';
    Vbold(iv).pinfo = [1,0,0];
    Vbold(iv).n = [iv 1];
end

Vbold = myspm_write_vol_4d(Vbold,bolddat);

funcdat = reshape(funcdat(:,:),s);
if contains(params.asl.splitaslbold,'meica'), funcdat = funcdat-bolddat+repmat(mean(funcdat,4),[1,1,1,s(4)]); end

Vasl = Vfunc;
for iv=1:numel(Vasl)
    Vasl(iv).fname = fullfile(ppparams.subperfdir,['f' ppparams.func(ie).prefix fname{1} '_asl.nii']);
    Vasl(iv).descrip = 'my_spmbatch - split asl';
    Vasl(iv).pinfo = [1,0,0];
    Vasl(iv).n = [iv 1];
end

Vasl = myspm_write_vol_4d(Vasl,funcdat);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_bold.nii'])};

ppparams.func(ie).funcfile = [fname{1} '_bold.nii'];
ppparams.func(ie).prefix = ['f' ppparams.func(ie).prefix];

ppparams.perf(ie).perffile = [fname{1} '_bold.nii'];
ppparams.perf(ie).prefix = ['f' ppparams.func(ie).prefix];

clear funcdat bolddat Vfunc Vbold Vasl
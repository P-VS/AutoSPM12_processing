function [ppparams,keepfiles,delfiles] = fmri_dune(ppparams,params,keepfiles,delfiles)

%% echte code
[dune_path,~,~] = fileparts(mfilename('fullpath'));
dune_path = fullfile(dune_path,'dune');

funcfile_e1 = [ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile];

%[~,cmdout] = system(sprintf(['cd ' dune_path ' && python3 dune.py ' ppparams.subfuncdir ' ' funcfile_e1 ' ' ...
%                                    ppparams.func(1).jsonfile ' ' num2str(numel(ppparams.echoes)) ' ' ppparams.der_file]),'-echo');

if contains(ppparams.func(1).funcfile,'_echo-')
    nfname = split(ppparams.func(1).funcfile,'_echo-');
    ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];
end

ppparams.func = ppparams.func(1);
ppparams.func(1).prefix = ['cd' ppparams.func(1).prefix];

ppparams.echoes = 1;
ppparams.meepi = false;

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(1).prefix ppparams.func(1).funcfile])};
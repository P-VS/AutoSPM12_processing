function [ppparams,keepfiles,delfiles] = my_spmbatch_dune(ppparams,params,keepfiles,delfiles)

if contains(ppparams.func(1).funcfile,'_echo-')
    nfname = split(ppparams.func(1).funcfile,'_echo-');
    ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];
end

ppparams.func = ppparams.func(1);
ppparams.func(1).prefix = ['cd' ppparams.func(1).prefix];

if params.func.isaslbold
    nfname = split(ppparams.perf(1).perffile,'_echo-');
    ppparams.perf(1).perffile = [nfname{1} '_asl.nii'];

    ppparams.perf(1).prefix = ['cd' ppparams.perf(1).prefix];

    oldfile = fullfile(ppparams.subfuncdir,[ppparams.perf(1).prefix ppparams.perf(1).funcfile]);
    newfile = fullfile(ppparams.subperfdir,[ppparams.perf(1).prefix ppparams.perf(1).funcfile]);

    movefile(oldfile,newfile);
end

ppparams.echoes = 1;
ppparams.meepi = false;

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(1).prefix ppparams.func(1).funcfile])};
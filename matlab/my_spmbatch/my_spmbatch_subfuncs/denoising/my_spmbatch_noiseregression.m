function [funcdat,ppparams,keepfiles] = my_spmbatch_noiseregression(funcdat,ppparams,params,keepfiles,regkind)

if ~params.preprocess_functional, params.denoise.before_normalization = false; end

if contains(regkind,'bold')
    if ~isfield(ppparams,'noiseregresssor')
        if params.denoise.do_ICA_AROMA && isfield(ppparams,'nboldica_file')
            confounds = load(ppparams.nboldica_file);
        else
            if params.denoise.do_mot_derivatives && isfield(ppparams,'der_file')
                confounds = load(ppparams.der_file);
            elseif isfield(ppparams,'rp_file')
                confounds = load(ppparams.rp_file);
            else
                confounds = [];
            end
        
            if params.denoise.do_aCompCor && isfield(ppparams,'acc_file')
                acc_confounds = load(ppparams.acc_file);
    
                confounds = cat(2,confounds,acc_confounds);
            end
        end
    else
        confounds = ppparams.noiseregresssor;
    end
elseif contains(regkind,'fasl') 
    if isfield(ppparams,'naslica_file')
        confounds = load(ppparams.naslica_file);
    else
        confounds = [];
    end
end

if params.denoise.do_bpfilter
    jsondat = fileread(ppparams.func(ppparams.echoes(1)).jsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    bpfilter = [tr params.denoise.bpfilter(1:2)];
else
    bpfilter = [];
end

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

[funcdat,~] = fmri_cleaning(funcdat(:,:),params.denoise.polort,bpfilter,confounds,[],'restoremean','on');

funcdat = reshape(funcdat(:,:),s);
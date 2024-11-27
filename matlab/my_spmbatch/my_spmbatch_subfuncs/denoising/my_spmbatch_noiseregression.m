function funcdat = my_spmbatch_noiseregression(funcdat,ppparams,params,regkind)

if ~params.preprocess_functional, params.denoise.before_normalization = false; end

mask = spm_read_vols(spm_vol(ppparams.fmask));

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

tmp=find(mask>0);
mfuncdat = funcdat(tmp,:);
mean_mfuncdat = mean(mfuncdat,2);
std_mfuncdat = std(mfuncdat,[],2);

mfuncdat = mfuncdat - repmat(mean_mfuncdat,[1,s(end)]);
mfuncdat(std_mfuncdat>0,:) = mfuncdat(std_mfuncdat>0,:) ./ repmat(std_mfuncdat(std_mfuncdat>0),[1,s(end)]);

if contains(regkind,'bold')
    if ~isfield(ppparams,'noiseregresssor')
        if params.denoise.do_ICA_AROMA
            if isfield(ppparams,'nboldica_file'), confounds = load(ppparams.nboldica_file); else confounds = []; end
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

[mfuncdat,~] = fmri_cleaning(mfuncdat(:,:),params.denoise.polort,[],confounds,[],'restoremean','on');

if params.denoise.do_bpfilter
    jsondat = fileread(ppparams.func(ppparams.echoes(1)).jsonfile);
    jsondat = jsondecode(jsondat);

    tr = jsondat.RepetitionTime;

    mfuncdat = filter_chunk(mfuncdat(:,:), params.denoise.bpfilter(1), params.denoise.bpfilter(2), tr);
end

mfuncdat = mfuncdat .* repmat(std_mfuncdat,[1,s(end)]);
mfuncdat = mfuncdat + repmat(mean_mfuncdat,[1,s(end)]);

funcdat(tmp,:) = mfuncdat;

funcdat = reshape(funcdat(:,:),s);

clear mfuncdat mean_mfuncdat std_mfuncdat tmp


function fil_calc = filter_chunk(vol_calc, filter_lower, filter_upper, filter_tr)

int_array = get_intervals(1, size(vol_calc, 2), 1024);
fil_calc = zeros(size(vol_calc), 'single');
for n = 1:size(int_array, 1)
    ind_vec = int_array(n, 1):int_array(n, 2);
    fil_calc(:, ind_vec) = brant_Filter_FFT_Butterworth(vol_calc(:, ind_vec), filter_lower, filter_upper, 1 / filter_tr, 0);
end

function int_array = get_intervals(start_pt, end_pt, len_int)
% works only for integer
% start_pt: start point
% end_pt: end point
% len_int: length of interval

assert(start_pt < end_pt);

num_pts = end_pt - start_pt + 1;

num_blk = ceil(double(num_pts) / double(len_int));

int_tmp = start_pt:len_int:end_pt;
int_array = zeros(num_blk, 2);
for m = 1:num_blk
    if m == num_blk
        int_array(m, :) = [int_tmp(m), min(end_pt, int_tmp(m) + len_int - 1)];
    else
        int_array(m, :) = [int_tmp(m), int_tmp(m) + len_int - 1];
    end
end
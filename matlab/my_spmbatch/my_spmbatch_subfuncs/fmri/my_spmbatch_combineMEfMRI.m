function [ppparams,delfiles] = my_spmbatch_combineMEfMRI(ppparams,params,delfiles)

fprintf('Combine echoes\n')

%% Loading the fMRI time series and deleting dummy scans

nechoes = numel(ppparams.echoes);

for ie=1:nechoes
    jsondat = fileread(ppparams.func(ppparams.echoes(ie)).jsonfile);
    jsondat = jsondecode(jsondat);

    te(ie) = 1000.0*jsondat.EchoTime;

    tefunc{ie}.Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ppparams.echoes(ie)).prefix ppparams.func(ppparams.echoes(ie)).funcfile]));
end

tdim = numel(tefunc{1}.Vfunc);
voldim = tefunc{1}.Vfunc(1).dim;

if params.func.isaslbold, params.func.combination='dyn_T2star'; end

spm_progress_bar('Init',tdim,'Combine TE images','volumes completed');

nvols = params.loadmaxvols;
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end
    tefuncdat = zeros(voldim(1),voldim(2),voldim(3),nvols,nechoes);
    
    for ie=1:nechoes
        tefuncdat(:,:,:,:,ie) = spm_read_vols(tefunc{ie}.Vfunc(ti:ti+nvols-1));
    end

    switch params.func.combination
        case 'average'
            funcdat = sum(tefuncdat,5) ./ nechoes;
    
        case 'TE_weighted'
        
            sum_weights = sum(te,'all');

            funcdat = zeros(voldim(1),voldim(2),voldim(3),nvols);
            
            for ie=1:nechoes
                functidat = tefuncdat(:,:,:,:,ie);
                functidat = functidat .* te(ie);
                functidat = functidat ./ sum_weights;
        
                funcdat = funcdat+functidat; 

                clear functidat
            end
    
        case 'T2star_weighted'
            if ti==1
                tefdat = mean(tefuncdat,4);

                t2star = make_t2star_map(tefdat,te);
                t2star = reshape(t2star,[voldim(1)*voldim(2)*voldim(3),1]);
                t2star_ind = find(t2star>0);

                clear tefdat

                weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);
            
                for ne=1:nechoes
                    weights(t2star_ind,ne) = repmat(-te(ne),numel(t2star_ind),1) ./ t2star(t2star_ind,1);
                    weights(t2star_ind,ne) = exp(weights(t2star_ind,ne));
                    weights(t2star_ind,ne) = repmat(te(ne),numel(t2star_ind),1) .* weights(t2star_ind,ne);
                end
            
                %weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes]);
                
                sum_weights = sum(weights,2);
                weights_mask = find(sum_weights>0);

                clear t2star_ind t2star
            end

            funcdat = zeros(voldim(1),voldim(2),voldim(3));

            for ne=1:nechoes
                functidat = reshape(tefuncdat(:,:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),nvols]);
                functidat = functidat .* repmat(weights(:,ne),[1 nvols]);
                functidat(weights_mask,:) = functidat(weights_mask,:) ./ repmat(sum_weights(weights_mask),[1 nvols]);
        
                funcdat = funcdat+reshape(functidat,[voldim(1),voldim(2),voldim(3),nvols]); 

                clear functidat
            end 

        case 'dyn_T2star'
                funcdat = make_t2star_map(tefuncdat,te);
                
    end

    Vout = tefunc{1}.Vfunc(ti:ti+nvols-1);
    rmfield(Vout,'pinfo');
    
    if ti==1, nfname = split(ppparams.func(1).funcfile,'_echo-'); end

    for iv=1:nvols
        Vout(iv).fname = fullfile(ppparams.subfuncdir,['c' ppparams.func(1).prefix nfname{1} '_bold.nii']);
        Vout(iv).descrip = 'my_spmbatch - combine echoes';
        Vout(iv).dt = [spm_type('float32'),spm_platform('bigend')];
        Vout(iv).n = [ti+iv-1 1];
        Vout(iv) = spm_write_vol(Vout(iv),funcdat(:,:,:,iv));
    end

    clear tefuncdat funcdat Vfunc 

    spm_progress_bar('Set',ti);
end    
spm_progress_bar('Clear');

nfname = split(ppparams.func(1).funcfile,'_echo-');

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['c' ppparams.func(1).prefix nfname{1} '_bold.nii']);};

ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];
ppparams.func(1).prefix = ['c' ppparams.func(1).prefix];

%% _______________________________________________________________________
function t2star = make_t2star_map(tefuncdat,te)
%based on https://github.com/jsheunis/fMRwhy/tree/master

voldim = size(tefuncdat);
if numel(voldim)<5, tefuncdat = reshape(tefuncdat,[voldim(1),voldim(2),voldim(3),1,voldim(4)]); end
voldim = size(tefuncdat);

nechoes = numel(te);

mask = zeros([voldim(1),voldim(2),voldim(3)]);

for ie=1:nechoes
    iemask = my_spmbatch_mask(tefuncdat(:,:,:,:,ie));
    mask = mask + iemask;

    clear iemask
end

mask_ind = find(mask>nechoes-0.5);

t2star = zeros(voldim(1),voldim(2),voldim(3),voldim(4));

for ti=1:voldim(4)
    % Create "design matrix" X
    X = horzcat(ones(nechoes,1), -te(:));
    
    tempt2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
    
    Y=[];
    for ne=1:nechoes
        temptefuncdat = reshape(tefuncdat(:,:,:,ti,ne),[voldim(1)*voldim(2)*voldim(3),1]);
        Y=[Y;reshape(temptefuncdat(mask_ind,1),[1,numel(mask_ind)])];
    end
    Y = max(Y, 1e-11);
    
    % Estimate "beta matrix" by solving set of linear equations
    beta_hat = pinv(X) * log(Y);
     % Calculate S0 and T2star from beta estimation
    T2star_fit = beta_hat(2, :); %is R2*
    
    tempt2star(mask_ind) = T2star_fit;
    tempt2star(mask_ind) = T2star_fit;
    zeromask = (tempt2star>0);
    
    tempt2star(zeromask) = 1 ./ tempt2star(zeromask);
    
    t2star_perc = prctile(tempt2star,99.5,'all');
    T2star_thresh_max = 10 * t2star_perc; %same as tedana
    tempt2star(tempt2star>T2star_thresh_max) = t2star_perc;
    t2star(:,:,:,ti) = reshape(tempt2star,[voldim(1),voldim(2),voldim(3)]);
    
    clear temptefuncdat tempt2star T2star_fit Y beta_hat
end

if numel(size(tefuncdat))<5, t2star = reshape(t2star,[voldim(1),voldim(2),voldim(3)]); end

clear mask mask_ind
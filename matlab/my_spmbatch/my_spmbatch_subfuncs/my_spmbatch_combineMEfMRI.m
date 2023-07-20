function [funcdat,Vout] = my_spmbatch_combineMEfMRI(tefuncdata,ppparams,params)
%% Loading the fMRI time series and deleting dummy scans

nechoes = numel(params.echoes);

for i=1:nechoes
    funcjsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring num2str(params.echoes(i)) '.json']);

    jsondat = fileread(funcjsonfile);
    jsondat = jsondecode(jsondat);

    te(i) = 1000.0*jsondat.EchoTime;

    Vfunc = tefuncdata{params.echoes(i)}.Vfunc;

    if i==1
        voldim = Vfunc.dim;
        tefuncdat = zeros(voldim(1),voldim(2),voldim(3),numel(Vfunc),nechoes);
    end

    tefuncdat(:,:,:,:,i) = tefuncdata{params.echoes(i)}.data;
end

funcdat = zeros(voldim(1),voldim(2),voldim(3),numel(Vfunc));

switch params.combination
    case 'average'
        funcdat = sum(tefuncdat,5) ./ nechoes;

    case'TE_weighted'
    
        sum_weights = sum(te,'all');
        
        for ti=1:numel(Vfunc)
            for ne=1:nechoes
                functidat = tefuncdat(:,:,:,ti,ne);
                functidat = functidat .* te(ne);
                functidat = functidat ./ sum_weights;
        
                funcdat(:,:,:,ti) = funcdat(:,:,:,ti)+functidat; 
            end
        end

    case'T2_weighted'
    
        mask = my_spmbatch_mask(tefuncdat(:,:,:,:,1));
        mask_ind = find(mask>0);

        %based on https://github.com/jsheunis/fMRwhy/tree/master
        for ti=1:numel(Vfunc)
    
            tifuncdat = reshape(tefuncdat(:,:,:,ti,:),[voldim(1),voldim(2),voldim(3),nechoes]);

            % Create "design matrix" X
            X = horzcat(ones(nechoes,1), -te(:));
        
            t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
        
            Y=[];
            for ne=1:nechoes
                temptefuncdat = reshape(tifuncdat(:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),1]);
                Y=[Y;reshape(temptefuncdat(mask_ind,1),[1,numel(mask_ind)])];
            end
            Y = max(Y, 1e-11);
        
            % Estimate "beta matrix" by solving set of linear equations
            beta_hat = pinv(X) * log(Y);
             % Calculate S0 and T2star from beta estimation
            T2star_fit = beta_hat(2, :); %is R2*
        
            T2star_thresh_min = 1/1500; % arbitrarily chosen, same as tedana
            I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
            T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
        
            t2star(mask_ind) = T2star_fit;
            
            weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);
        
            for ne=1:nechoes
                weights(:,ne) = repmat(-te(ne),voldim(1)*voldim(2)*voldim(3),1) .* t2star(:,1);
                weights(:,ne) = exp(weights(:,ne));
                weights(:,ne) = repmat(te(ne),voldim(1)*voldim(2)*voldim(3),1) .* weights(:,ne);
            end
        
            weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes]);
            
            sum_weights = sum(weights,4);
            weights_mask = find(sum_weights>0);

            for ne=1:nechoes
                functidat = tefuncdat(:,:,:,ti,ne);
                functidat = functidat .* weights(:,:,:,ne);
                functidat(weights_mask) = functidat(weights_mask) ./ sum_weights(weights_mask);
        
                funcdat(:,:,:,ti) = funcdat(:,:,:,ti)+functidat; 
            end 
        end

    case 'dyn_T2star'

        mask = my_spmbatch_mask(tefuncdat(:,:,:,:,1));
        mask_ind = find(mask>0);

        %based on https://github.com/jsheunis/fMRwhy/tree/master
        for ti=1:numel(Vfunc)
    
            tifuncdat = reshape(tefuncdat(:,:,:,ti,:),[voldim(1),voldim(2),voldim(3),nechoes]);

            % Create "design matrix" X
            X = horzcat(ones(nechoes,1), -te(:));
        
            t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
        
            Y=[];
            for ne=1:nechoes
                temptefuncdat = reshape(tifuncdat(:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),1]);
                Y=[Y;reshape(temptefuncdat(mask_ind,1),[1,numel(mask_ind)])];
            end
            Y = max(Y, 1e-11);
        
            % Estimate "beta matrix" by solving set of linear equations
            beta_hat = pinv(X) * log(Y);
             % Calculate S0 and T2star from beta estimation
            T2star_fit = beta_hat(2, :); %is R2*
        
            T2star_thresh_min = 1/1500; % arbitrarily chosen, same as tedana
            I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
            T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
        
            t2star(mask_ind) = T2star_fit;
            zeromask = (t2star>0);

            t2star(zeromask) = 1 ./ t2star(zeromask);

            funcdat(:,:,:,ti) = reshape(t2star,[voldim(1),voldim(2),voldim(3)]);
        end

end

Vfunc = tefuncdata{ppparams.echoes(1)}.Vfunc;

[fpath,fname,~] = fileparts(Vfunc(1).fname);
nfname = split(fname,'bold_e');

Vout = Vfunc;

for j=1:numel(Vout)
    Vout(j).fname = fullfile(fpath,['c' ppparams.prefix nfname{1} 'bold.nii']);
    Vout(j).descrip = 'my_spmbatch - combine echoes';
    Vout(j).pinfo = [1,0,0];
    Vout(j).dt = [spm_type('float32'),spm_platform('bigend')];
    Vout(j).n = [j 1];
end

Vout = myspm_write_vol_4d(Vout,funcdat);
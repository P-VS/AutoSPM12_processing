function [Vfunc,funcdat,funcfile] = my_spmbatch_readMEfMRI(subfuncstring,subfmridir,numdummy,params,readvols)
%% Loading the fMRI time series and deleting dummy scans
fprintf('Reading the data \n')

for i=1:params.nechoes
    funcjsonfile = fullfile(subfmridir,[subfuncstring num2str(i) '.json']);
    funcfile = fullfile(subfmridir,[subfuncstring num2str(i) '.nii']);

    jsondat = fileread(funcjsonfile);
    jsondat = jsondecode(jsondat);

    te(i) = 1000.0*jsondat.EchoTime;

    Vfunc = spm_vol(funcfile);

    if params.reorient
        transfile = fullfile(subfmridir,[subfuncstring '1_reorient.mat']);
        if isfile(transfile)
            load(transfile,'M')
            transM = M;
        else
            transM = eye(4);
        end

        MM = Vfunc(1).private.mat0;

        Vfunc = my_reset_orientation(Vfunc,transM * MM);
    end
    
    Vfunc = Vfunc(numdummy+1:end);
    if readvols<Inf
        Vfunc = Vfunc(1:readvols); %Only for a quick test of the batch script
    end

    efuncdat = spm_read_vols(Vfunc);

    if i==1
        voldim = Vfunc.dim;
        tefuncdat = zeros(voldim(1)*voldim(2)*voldim(3),numel(Vfunc),params.nechoes);
    end

    tefuncdat(:,:,i) = reshape(efuncdat,[voldim(1)*voldim(2)*voldim(3),numel(Vfunc)]);
end

weights = zeros(voldim(1)*voldim(2)*voldim(3),params.nechoes);

if params.combination=='TE_weighted'
    for ne=1:params.nechoes
        weights(:,ne) = repmat(te(ne),voldim(1)*voldim(2)*voldim(3),1);
    end
elseif params.combination=='T2_weighted'
    %based on https://github.com/jsheunis/fMRwhy/tree/master

    % Create "design matrix" X
    X = horzcat(ones(params.nechoes,1), -te(:));

    t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);

    thres = 0.02*max(tefuncdat(:,1,1),[],'all');
    mask_ind = find(tefuncdat(:,1,1)>thres);

    Y=[];
    for ne=1:params.nechoes
        Y=[Y;reshape(tefuncdat(mask_ind,1,ne),[1,numel(mask_ind)])];
    end
    Y = max(Y, 1e-11);

    % Estimate "beta matrix" by solving set of linear equations
    beta_hat = pinv(X) * log(Y);
     % Calculate S0 and T2star from beta estimation
    T2star_fit = 1 ./ beta_hat(2, :);

    T2star_thresh_max = 5000; % arbitrarily chosen
    T2star_thresh_min = 0; % arbitrarily chosen, same as tedana
    I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
    T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
    I_T2star_max =  T2star_fit(:) >= T2star_thresh_max; % vector of voxels in mask where T2star value is higher than threshold_max
    T2star_fit(I_T2star_max) = 0; % set to zero, different from tedana where they set it to threshold_max value

    t2star(mask_ind)=T2star_fit;
    wmask = find(t2star>0);

    for ne=1:params.nechoes
        weights(wmask,ne)=te(ne).*exp(-te(ne)./t2star(wmask,1));
    end
end

sum_weights = sum(weights,2);
weights_mask = find(sum_weights>0);

funcdat = zeros(voldim(1)*voldim(2)*voldim(3),numel(Vfunc));

for ti=1:numel(Vfunc)
    for ne=1:params.nechoes
        funcdat(weights_mask,ti) = funcdat(weights_mask,ti)+(reshape(tefuncdat(weights_mask,ti,ne),[numel(weights_mask),1]) .* weights(weights_mask,ne));
    end

    funcdat(weights_mask,ti) = funcdat(weights_mask,ti) ./ sum_weights(weights_mask);
end

funcdat = reshape(funcdat,[voldim(1),voldim(2),voldim(3),numel(Vfunc)]);

fname=Vfunc(1).fname;
nfname = split(fname,['_e' num2str(params.nechoes)]);

for j=1:numel(Vfunc)
    Vfunc(j).fname = [nfname{1} '.nii'];
    Vfunc(j).pinfo = [];
end

funcfile = [nfname{1} '.nii'];
function [ppparams,delfiles,keepfiles] = my_spmbatch_fasl_cbfmapping(ppparams,params,delfiles,keepfiles)

%used CBF model ass in Alsop et al. Recommended Implementation of Arterial Spin-Labeled Perfusion MRI for Clinical Applications: 
% %A Consensus of the ISMRM Perfusion Study Group and the European
% Consortium for ASL in Dementia. Magnetic Resonance in Medicine 73:102–116 (2015)
%(https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.25197?src=getftr)

GM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile));
WM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile));
CSF = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile));

gmim = spm_read_vols(GM);
wmim = spm_read_vols(WM);
csfim = spm_read_vols(CSF);

Vdeltam=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).deltamprefix ppparams.perf(1).deltamfile]));

nfname = split(ppparams.perf(1).deltamfile,'_deltam');

tdim = numel(Vdeltam);
voldim = Vdeltam(1).dim;

jsondat = fileread(ppparams.func(1).jsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;
te = 1000*jsondat.EchoTime;
if isfield(jsondat,'LabelingDuration'), LD = jsondat.LabelingDuration; else LD = params.asl.LabelingDuration; end
if isfield(jsondat,'PostLabelDelay'), PLD = jsondat.PostLabelDelay; else PLD = params.asl.PostLabelDelay; end

if isfield(jsondat,'SliceTiming'), SliceTimes = jsondat.SliceTiming; else SliceTimes = []; end
if ~(numel(SliceTimes)==voldim(3)), SliceTimes = []; end

if isempty(SliceTimes)
    if isfield(jsondat,'MultibandAccelerationFactor')
        hbf = jsondat.MultibandAccelerationFactor; 
        nslex = ceil(voldim(3)/hbf);
        isl = zeros([1,nslex]);
        isl(1:2:nslex)=[0:1:(nslex-1)/2];
        isl(2:2:nslex)=[ceil(nslex/2):1:nslex-1];
        isl=repmat(isl,[1,hbf]);
    else 
        isl = [1:2:voldim(3) 2:2:voldim(3)];
        nslex = voldim(3);
    end

    TA = tr-LD-PLD;
    SliceTimes = isl*TA/(nslex-1);
else
    TA = tr-LD-PLD;

    SliceTimes = SliceTimes * TA/tr;
end

SlicePLD = PLD+SliceTimes;

T1a = 1.650; %longitudinal relaxation time of arterial blood
lambda = 0.9; %blood-brain partition coefficient for gray matter
alpha = 0.85; %laeling efficiency

vol_PLD = reshape(repmat(SlicePLD,[voldim(1)*voldim(2),1]),voldim);

Vm0=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile]));
m0vol = spm_read_vols(Vm0);

%ccorrect M0 for T1 effects
% The T1 values used, are the averaged T1 values reported in the review of 
% Bojorquez et al. 2017. What are normal relaxation times of tissues at 3 T? Magnetic Resonance Imaging 35:69-80
% (https://mri-q.com/uploads/3/4/5/7/34572113/normal_relaxation_times_at_3t.pdf)

T1gm = 1.459;
T1wm = 0.974;
T1csf = 4.190;

T1dat = (T1gm * gmim + T1wm * wmim + T1csf * csfim);

corr_T1 = zeros(voldim);
corr_T1(T1dat>0) = 1 ./ (1-exp(-tr./T1dat(T1dat>0)));
m0vol = m0vol .* corr_T1;

mask = my_spmbatch_mask(m0vol);

cm0vol = 2*alpha*m0vol*T1a.*(exp(-vol_PLD/T1a)-exp(-(LD+vol_PLD)/T1a));
cm0vol = reshape(cm0vol,[voldim(1)*voldim(2)*voldim(3),1]);

clear gmim wmim csfim GM WM CSF

nvols = params.loadmaxvols;
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end

    fprintf(['CBF vols: ' num2str(ti) '-' num2str(ti+nvols-1) '\n'])
    
    cbfdata = spm_read_vols(Vdeltam(ti:ti+nvols-1));
    cbfdata = reshape(cbfdata,[voldim(1)*voldim(2)*voldim(3),nvols]);

    cbfdata = lambda*6000*cbfdata;
    cbfdata(cm0vol>0,:) = cbfdata(cm0vol>0,:)./repmat(cm0vol(cm0vol>0),[1,nvols]);
    cbfdata = reshape(cbfdata,[voldim(1),voldim(2),voldim(3),nvols]);

    cbfdata = cbfdata.*repmat(mask,[1,1,1,nvols]);
    cbfdata(cbfdata<-40) = 0;
    cbfdata(cbfdata>150) = 0;

    Vout = Vdeltam(1);
    rmfield(Vout,'pinfo');
    for iv=1:nvols
        Vout.fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).deltamprefix nfname{1} '_cbf.nii']);
        Vout.descrip = 'my_spmbatch - cbf';
        Vout.dt = [spm_type('float32'),spm_platform('bigend')];
        Vout.n = [ti+iv-1 1];
        Vout = spm_write_vol(Vout,cbfdata(:,:,:,iv));
    end

    clear cbfdata Vout
end

clear m0vol cm0vol vol_PLD mask

ppparams.perf(1).cbfprefix = ppparams.perf(1).deltamprefix;
ppparams.perf(1).cbffile = [nfname{1} '_cbf.nii'];
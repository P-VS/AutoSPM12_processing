function p_out = acid_realign(PS, PT, PO)

% =========================================================================
% The script allows the user to realign images in a slice-wise manner.
%
% Inputs:
%   PS    - source image
%   PT    - target image
%   PO    - other images
%   p_out - output directory
%
% Outputs:
%   p_out - output foldername
%   fname - output filename
%
% =========================================================================

interp_def = -4;
slice = 1;
subj  = 1;
cax   = [0 1];

% load in source image
VS = spm_vol(PS);

% load in target image
VT  = spm_vol(PT);
dm  = VT.dim;
vol = readvol(dm,VT,VT,interp_def,[],slice);

% load in other image(s)
if size(PO,1)
    VO = spm_vol(PO);
else
    VO = acid_load_4Dimage(PO);
end

% define output directory
keyword = 'REALIGN';
[path,fname,~] = spm_fileparts(VS.fname);

p_out = acid_bids(path,fname,keyword,1);

    
if ~isempty(VO)
    % save json file
    acid_save_json(VO(1), p_out, keyword);
    
    % save bvals and bvecs files
    acid_save_bvals_bvecs(VO(1), p_out, keyword);
end

% params contains the degree of freedom (max 12) for each slice.
% Scaling in each dimension is allowed for all slices.
params = zeros(dm(3),12);
params(:,7:9) = ones(dm(3),3);

figure(2); clf;
subplot(1,2,1)
imagesc(vol);
caxis(cax);
colormap gray
axis off
title('Reference')

alpha = [1e2 1e-2];
alpha2 = [1 1e-2];
a1Range = [-100 100];
alphaSTART1 = a1Range(1);
a2Range = [-0.5 0.5];

%% GUI OBJECTS
% sliders for Y (Ya1Slider - translation, Ya2Slider - scaling)
Ya1Slider = uicontrol('Style','slider','Min',a1Range(1),'Max',a1Range(2),...
    'SliderStep',[.25 .25]./(a1Range(2)-a1Range(1)),'Value',log10(alpha(1)),...
    'Position',[180 45 300 25],'Callback',@Callback_Ya1Slider);

Ya2Slider = uicontrol('Style','slider','Min',a2Range(1),'Max',a2Range(2),...
    'SliderStep',1e-1*[.25 .25]./(a2Range(2)-a2Range(1)),'Value',log10(1+alpha2(1)),...
    'Position',[180 15 300 25],'Callback',@Callback_Ya2Slider);

% sliders for X (Xa1Slider - translation, Xa2Slider - scaling)
Xa1Slider = uicontrol('Style','slider','Min',a1Range(1),'Max',a1Range(2),...
    'SliderStep',[.25 .25]./(a1Range(2)-a1Range(1)),'Value',log10(alpha(1)),...
    'Position',[720 45 300 25],'Callback',@Callback_Xa1Slider);

Xa2Slider = uicontrol('Style','slider','Min',a2Range(1),'Max',a2Range(2),...
    'SliderStep',1e-1*[.25 .25]./(a2Range(2)-a2Range(1)),'Value',log10(1+alpha2(1)),...
    'Position',[720 15 300 25],'Callback',@Callback_Xa2Slider);

% slider for slice selection
sliceSlider = uicontrol('Style','slider','Min',1,'Max',dm(3),...
    'SliderStep',[1 1]/(dm(3)-1),'Value',slice,...
    'Position',[80 270 40 400]);
set(sliceSlider,'Callback',@Callback_sliceSlider);

% button for applying transformation
cax1Button = uicontrol('Style','push',...
    'Position', [1100 20 130 60], 'FontWeight', 'bold','string','Apply Transformation','Callback',@Callback_applytransf);

% texts on the interface
sliceText = uicontrol('Style','text',...
    'String',['slice=' num2str(slice)],...
    'Position',[50 680 100 20], 'FontSize', 15, 'FontWeight', 'bold');
Ya1Text = uicontrol('Style','text',...
    'String',['y-position=' num2str(alpha(1)+alphaSTART1)],...
    'Position',[40 47 120 20], 'FontSize', 15, 'FontWeight', 'bold');
Ya2Text = uicontrol('Style','text',...
    'String',['y-scaling=' num2str(1+alpha2(1))],...
    'Position',[40 17 120 20], 'FontSize', 15, 'FontWeight', 'bold');
Xa1Text = uicontrol('Style','text',...
    'String',['x-position=' num2str(alpha(1)+alphaSTART1)],...
    'Position',[580 47 120 20], 'FontSize', 15, 'FontWeight', 'bold');
Xa2Text = uicontrol('Style','text',...
    'String',['x-scaling=' num2str(1+alpha2(1))],...
    'Position',[580 17 120 20], 'FontSize', 15, 'FontWeight', 'bold');

% plot source and target images
showReg();

%% FUNCTION DEFINITIONS
% function for plotting images
function showReg()
    subplot(1,2,1);
    title(sprintf('Target'), 'FontSize', 20, 'FontWeight', 'bold');
    hold on;
    vol1 = readvol(dm, VT, VT, interp_def, [], slice);
    imagesc(vol1);
    caxis(cax);
        
    subplot(1,2,2)
    imagesc(vol);
    title(sprintf('Source'), 'FontSize', 20, 'FontWeight', 'bold');
    caxis(cax);
    axis off
    subplot(1,2,1);
    hold on;
    contour(vol,'r');
end

% callback function for updating slice ()
function Callback_sliceSlider(hObject,eventdata)
    slicetmp = slice;
    slice = round(get(hObject,'Value'));
    deActivateSliders();
    paramstmp = params(slicetmp,:);
    params(slice,:) = paramstmp;
    vol2D = readvol(dm,VT,VT,interp_def,paramstmp,slice);
        
    subplot(1,2,1)
    imagesc(vol2D)
    colormap gray
    axis off
    title('Reference')
        
    vol = readvol(dm,VT,VS(subj),interp_def,paramstmp,slice);
    showReg()
    activateSliders();
    clear slicetmp paramstmp;
end

% callback function for translation in Y
function Callback_Ya1Slider(hObject,eventdata)
    alpha(1) = (get(hObject,'Value'));
    params(slice,2) = alpha(1);
    xpar = params(slice,:);
    deActivateSliders();
    vol = readvol(dm,VT,VS(subj),interp_def,xpar,slice);
    clc
    showReg()
    activateSliders();
end

% callback function for scaling in Y
function Callback_Ya2Slider(hObject,eventdata)
    alpha2(1) = (get(hObject,'Value'));
    params(slice,8) = 1+alpha2(1);
    xpar = params(slice,:);
    deActivateSliders();
    vol = readvol(dm,VT,VS(subj),interp_def,xpar,slice);
    showReg()
    activateSliders();
end

% callback function for translation in X
function Callback_Xa1Slider(hObject,eventdata)
    alpha(1) = (get(hObject,'Value'));
    params(slice,1) = alpha(1);
    xpar = params(slice,:);
    deActivateSliders();
    vol = readvol(dm,VT,VS(subj),interp_def,xpar,slice);
    showReg()
    activateSliders();
end

% callback function for scaling in X
function Callback_Xa2Slider(hObject,eventdata)
    alpha2(1) = (get(hObject,'Value'));
    params(slice,7) = 1+alpha2(1);
    xpar = params(slice,:);
    deActivateSliders();
    vol = readvol(dm,VT,VS(subj),interp_def,xpar,slice);
    showReg()
    activateSliders();
end

function activateSliders()
    set(Ya1Slider,'Enable','on')
    set(Ya2Slider,'Enable','on')
    set(Xa1Slider,'Enable','on')
    set(Xa2Slider,'Enable','on')
    set(sliceSlider,'Enable','on')
end

function deActivateSliders()
    set(Ya1Text,'String',['y-pos=' num2str(alpha(1))])
    set(Ya2Text,'String',['y-scale=' num2str(alpha2(1))])
    set(Xa1Text,'String',['x-pos=' num2str(alpha(1))])
    set(Xa2Text,'String',['x-scale=' num2str(alpha2(1))])
    set(sliceText,'String',['slice=' num2str(slice)])
    set(Ya1Slider,'Enable','off')
    set(Ya2Slider,'Enable','off')
    set(Xa1Slider,'Enable','off')
    set(Xa2Slider,'Enable','off')
    set(sliceSlider,'Enable','off')
    pause(.001);
end

function vol = readvol(dm,VT,VF,interp_def,params,zpos)
    if exist('params','var') && ~isempty(params)
        iM = inv(acid_spm_matrix(params));
        if(exist('zpos','var'))
            M = spm_matrix([0 0 -zpos 0 0 0 1 1 1]);
            M1 = inv(M*inv(VT.mat)*iM*VF.mat);
            vol = spm_slice_vol(VF,M1,dm(1:2),interp_def);
        else
            vol = zeros(dm);
            for zpos=1:dm(3)
                M = spm_matrix([0 0 -zpos 0 0 0 1 1 1]);
                M1 = inv(M*inv(VT.mat)*iM*VF.mat);
                vol(:,:,zpos) = spm_slice_vol(VF,M1,dm(1:2),interp_def);
            end
        end
    else
        if exist('zpos','var')
            M   = spm_matrix([0 0 -zpos 0 0 0 1 1 1]);
            M1  = inv(M*inv(VT.mat)*VF.mat);
            vol = spm_slice_vol(VF,M1,dm(1:2),interp_def);
        else
            vol = zeros(dm);
            for zpos=1:dm(3)
                M = spm_matrix([0 0 -zpos 0 0 0 1 1 1]);
                M1 = inv(M*inv(VT.mat)*VF.mat);
                vol(:,:,zpos) = spm_slice_vol(VF,M1,dm(1:2),interp_def);
            end
        end
    end
    if any(vol(:)>0)
        [y,x] = hist(vol(find(vol>0)),size(vol,1));
        cy    = cumsum(y);
        sz    = size(vol(find(vol>0)),1);
        
        tmp = find(cy<=sz*0.9);
        if(~isempty(tmp))
            THR = x(max(tmp));
        else
            THR = max(vol(:));
        end
    else
        THR = 1;
    end
    vol = vol/THR;
end

% callback function for applying transformation
function Callback_applytransf(hObject, eventdata)

    % apply transformation to source
    I_rsource = zeros(dm);
    for p = 1:dm(3)
        xpar = params(p,:);
        iM = inv(acid_spm_matrix(xpar));
        M = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(VT.mat)*iM*VS(subj).mat);
        tmp = spm_slice_vol(VS(subj),M,dm(1:2),interp_def);
        I_rsource(:,:,p) = tmp;
    end
        
    Vout = VT;
    Vout.fname = VS(subj).fname;

    % write realigned source image
    keyword = 'REALIGN';
    acid_write_vol(I_rsource, Vout, p_out, keyword, 'same');
        
    % apply transformation to other images
    sz = size(VO,1);
    for i = 1:sz
        
        I_rothers = zeros(dm);
        Vout(i,1) = VT;
        [p,f,e] = fileparts(VO(i).fname);
        Vout(i,1).fname = [p filesep f e];
        Vout(i,1).n = VO(i).n;
        Vout(i,1).descrip = VO(i).descrip;
        
        for p = 1:dm(3)
            xpar = params(p,:);
            iM   = inv(acid_spm_matrix(xpar));
            M    = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(VT.mat)*iM*VO(i).mat);
            tmp  = spm_slice_vol(VO(i),M,dm(1:2),interp_def);
            I_rothers(:,:,p) = tmp;
        end
            
        I_out(:,:,:,i) = I_rothers;
        
    end
    
    % write realigned images
    if size(VO,1)==1
        V_out = acid_write_vol(I_out, Vout, p_out, keyword, 'same');
    elseif strcmp(VO(1).fname,VO(2).fname)
        V_out = acid_write_vol(I_out, Vout(1), p_out, keyword, 'same');
    else
        for i=1:size(VO,1)
            V_out(i) = acid_write_vol(I_out(:,:,:,i), Vout(i), p_out, keyword, 'same');
        end
    end
    
    % save mat file
    fname_params = acid_bids_filename(VS,[keyword '-params'],'','.mat');
    fname_params = [p_out filesep fname_params];
    save(fname_params, 'params');
    
end
end
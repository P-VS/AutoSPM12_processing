function acid_dwi_series_browser(P, bvals, bvecs)

% =========================================================================
% The scripts is designed to browse through a dMRI dataset.
%
% Inputs:
%   P     - dMRI dataset
%   bvals - b values
%   bvecs - b vectors
%   p_out - output directory
% =========================================================================

slice = 1;
vol   = 1;
k     = 1;
idx   = [];

% check input
if size(P,1)==0
    P = char(cfg_getfile(Inf,'IMAGE','Select images','','','.*'));
end

% load in 4D dataset
V = acid_load_4Dimage(P);
I = spm_read_vols(V);

% check bvals/bvecs input
if ~isempty(bvals)
    dummy_bvals = 1;
    if size(bvals,2)~=size(V,1), error('The number of source images and b-value entries has to be the same.'); end
else
    dummy_bvals = 0;
end
if ~isempty(bvecs)
    dummy_bvecs = 1;
    if size(bvecs,2)~=size(V,1), error('The number of source images and b-vectors entries has to be the same.'); end
else
    dummy_bvecs = 0;
end
    
% get output directory
[p_out,~,~] = fileparts(P(1,:));


% get output filename
keyword = 'labels';
fname = acid_bids_filename(V, keyword, '', '.mat');
fname = [p_out filesep fname];

% determine axis length ratios (for plotting)
dm = V(1).dim;
if dm(1)>dm(2)
    ratio = [round(dm(1)/dm(2)),1,1];
else
    ratio = [1,round(dm(2)/dm(1)),1];
end

figure(2)
PlotSlice();

%% GUI OBJECTS

slider_vols = uicontrol('Style','slider','Min',1,'Max',size(V,1),...
    'SliderStep',[1/(size(V,1)-1) 1/(size(V,1)-1)],'Value',vol,...
    'Position',[40 400 40 400],'Callback',@Callback_slider_vols);

slider_slices = uicontrol('Style','slider','Min',1,'Max',size(I,3),...
    'SliderStep',[1/(size(I,3)-1) 1/(size(I,3)-1)],'Value',slice,...
    'Position',[140 400 40 400],'Callback',@Callback_slider_slices);

label_vols = uicontrol('Style','text','String',['Volume = ' num2str(vol)],...
    'Position',[15 870 90 15], 'FontSize', 10, 'FontWeight', 'bold');

label_slices = uicontrol('Style','text','String',['Slice = ' num2str(slice)],...
    'Position',[115 870 90 15], 'FontSize', 10, 'FontWeight', 'bold');

button_outlier = uicontrol('Style','push','Position', [35 300 150 80],...
    'string','Label', 'FontSize', 15, 'FontWeight', 'bold','Callback',@Callback_button_outlier);

%% FUNCTION DEFINITIONS

% function for plotting slice
function PlotSlice()
    imagesc(imrotate(I(:,:,slice,vol),90));
    colormap('gray')
    if (dummy_bvals && dummy_bvecs)
        title(['dMRI image: bval=' num2str(bvals(vol)) '  bvec=[' num2str(bvecs(1,vol)) ', ' num2str(bvecs(2,vol)) ', ' num2str(bvecs(3,vol)) ']'], 'FontSize', 18, 'FontWeight', 'bold')
    elseif dummy_bvals
        title(['dMRI image: bval=' num2str(bvals(vol))], 'FontSize', 18, 'FontWeight', 'bold')
    elseif dummy_bvecs
        title(['dMRI image: bvec=[' num2str(bvecs(1,vol)) ', ' num2str(bvecs(2,vol)) ', ' num2str(bvecs(3,vol)) ']'], 'FontSize', 18, 'FontWeight', 'bold')
    else
        title('dMRI image')
    end
    % axis square
    set(gca, 'fontsize', 15, 'TickDir', 'out', 'TitleFontSizeMultiplier', 1.5)
    %xlabel('FontSize', 18)
    %ylabel('FontSize', 18)
    pbaspect(ratio)
 end

% callback function for updating volume
function Callback_slider_vols(hObject,eventdata,handles)
    vol = round(get(hObject,'Value'));
    deActivateSliders();
    PlotSlice();
    ActivateSliders()
end

% callback function for updating slice
function Callback_slider_slices(hObject,eventdata,handles)
    slice = round(get(hObject,'Value'));
    deActivateSliders();
    PlotSlice();
    ActivateSliders()
end

function deActivateSliders()
    set(label_vols,'String',['Vol = ' num2str(vol)])
    set(label_slices,'String',['Slc = ' num2str(slice)])
    set(slider_vols,'Enable','off')
    set(slider_slices,'Enable','off')
    pause(.001);
end

function ActivateSliders()
    set(slider_vols,'Enable','on')
    set(slider_slices,'Enable','on')
    pause(.001);
end

% callback function for labeling slices
function Callback_button_outlier(hObject,eventdata,handles)
    
    if k==1
        idx(:,k) = [vol,slice];
        save(fname,'idx')
        k = k + 1;
    elseif ~ismember([vol,slice],idx','rows')
        idx(:,k) = [vol,slice];
        idx = sortrows(idx')';
        save(fname,'idx')
        k = k + 1;
    end
end

end
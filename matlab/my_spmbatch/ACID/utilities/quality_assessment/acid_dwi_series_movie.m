function acid_dwi_series_movie(P_ref, P1, P2, P3, dummy_allslices, dummy_order, slice, intv_ref, intv_src, dummy_movie, Pause,dummy_contour)

% determine number of datasets
if isempty(P1)
    error('Select at least one dataset')
elseif isempty(P2)
    N = 1;
elseif isempty(P3)
    N = 2;
else
    N = 3;
end


% Number of lines for contour plot
N_lines = 3;

% load in reference image
V_ref = spm_vol(P_ref);
dm = V_ref.dim;
I_ref = spm_read_vols(V_ref);


% check for correctness of slice

if isempty(slice)
    slice = ceil(dm(3)/2);
end

if (dummy_allslices && (slice<1 || slice>dm(3)))
    error('The specified slice number does lies outside the image');
end

% load in 1. dataset
V1 = acid_load_4Dimage(P1);
x = V1.dim;
slice_for_ref = round(slice/x(3)*dm(3));
I1 = spm_read_vols(V1);

% load in 2. dataset
if N > 1
    if (size(P2,1)~=size(P1,1))
        error('The size of the first and second dataset is not the same.')
    end

    V2 = acid_load_4Dimage(P2);
    I2 = spm_read_vols(V2);
    if V2(1).dim(3) ~= V1(1).dim(3)
        error('The number of slices in the first and second dataset is not the same.');
    end
end

% load in 3. dataset
if N == 3
    if (size(P3,1)~=size(P1,1))
        error('The size of the first and third dataset is not the same.');
    elseif (size(P3,1)~=size(P2,1))
        error('The size of the second and third dataset is not the same.');
    end
    V3 = acid_load_4Dimage(P3);
    I3 = spm_read_vols(V3);

    if V1(1).dim(3) ~= V3(1).dim(3)
        error('The number of slices in the first and third dataset is not the same.');
    elseif V2(1).dim(3) ~= V3(1).dim(3)
        error('The number of slices in the second and third dataset is not the same.');
    end
end

if dummy_contour == 0

    I_contour = I_ref;
    contour_text = 'Contours: Reference image ';

elseif dummy_contour == 1
    contour_text = 'Contours: Dataset 1';

    I_contour = I1;

elseif dummy_contour == 2
    contour_text = 'Contours: Dataset 2';

    I_contour = I2;

elseif dummy_contour == 3
    contour_text = 'Contours: Dataset 3';

    I_contour = I3;

elseif dummy_contour == 4
    contour_text = '';
    I_contour = I_ref;
    N_lines = 0;

end



fig1 = figure;
set(fig1,'units','normalized','outerposition', [0 0 1 1]);
colormap gray;
set(fig1, 'Name', 'QC of motion correction');
if dummy_movie
 
    % get output filename
    keyword = 'MOVIE';
    switch N
        case 1
            keyword_file = 'ONE-DATASET-MOVIE';
            fname_ref = acid_bids_filename(V_ref, keyword_file, '', '.avi');
            [path,~,~] = spm_fileparts(V_ref.fname);
            % fname = [fname_ref '--' fname_1];
        case 2
            keyword_file = 'TWO-DATASETS-MOVIE';
            fname_ref = acid_bids_filename(V_ref, keyword_file, '', '.avi');
            [path,~,~] = spm_fileparts(V_ref.fname);
            % fname = [fname_ref '--' fname_1 '--' fname_2];
        case 3
            keyword_file = 'THREE-DATASETS-MOVIE';
            fname_ref = acid_bids_filename(V_ref, keyword_file, '', '.avi');            
            [path,~,~] = spm_fileparts(V_ref.fname);
            % fname = [fname_ref '--' fname_1 '--' fname_2 '--' fname_3];
    end
    
    % % get output directory
    % 
    % switch N
    %     case 1
    %         [p_out,~,~] = fileparts(P1(1,:));
    %     case 2
    %         [p_out,~,~] = fileparts(P2(1,:));
    %     case 3
    %         [p_out,~,~] = fileparts(P3(1,:));
    % end
    % [path,fname,~] = spm_fileparts(V_ref.fname);
    p_out = acid_bids(path,fname_ref,keyword,1);

    
    fname = [p_out filesep fname_ref];
  
    % starting a movie object   
    cd(p_out)
    aviobj = VideoWriter(fname);
    open(aviobj);
end


switch N

    % 1 dataset
    case 1
        
        % all slices
        if dummy_allslices
            if ~dummy_order
                for k=1:dm(3)
                    for i=1:size(I1,4)                        
                        subplot(1,2,1);imagesc(rot90(I_ref(:,:,k)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',k)}, 'FontSize', 18, 'FontWeight', 'bold'); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,2,2);imagesc(rot90(I1(:,:,k,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',k,i)}, 'FontSize', 18, 'FontWeight', 'bold'); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        
                        sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                        
                        pause(Pause)
                        if dummy_movie
                            frame = getframe(gcf);
                            writeVideo(aviobj,frame);
                        end
                    end
                end
            else
                for i=1:size(I1,4)
                    for k=1:dm(3)
                        subplot(1,2,1);imagesc(rot90(I_ref(:,:,k)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',k)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,2,2);imagesc(rot90(I1(:,:,k,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                        pause(Pause)
                        if dummy_movie
                            frame = getframe(gcf);
                            writeVideo(aviobj,frame);
                        end
                    end
                end
            end

        % single slice
        else
            for i=1:size(I1,4)
                subplot(1,2,1);imagesc(rot90(I_ref(:,:,slice)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',slice)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                subplot(1,2,2);imagesc(rot90(I1(:,:,slice,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',slice,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                pause(Pause)
                if dummy_movie
                    frame = getframe(gcf);
                    writeVideo(aviobj,frame);
                end
            end
        end

    % 2 datasets
    case 2
        
        % all slices
        if dummy_allslices
            if ~dummy_order
                for k=1:dm(3)
                    for i=1:size(I1,4)
                        subplot(1,3,1);imagesc(rot90(I_ref(:,:,k)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',k)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,3,2);imagesc(rot90(I1(:,:,k,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,3,3);imagesc(rot90(I2(:,:,k,i)),intv_src); axis off; title({'Images 2', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                        pause(Pause)
                        if dummy_movie
                            frame = getframe(gcf);
                            writeVideo(aviobj,frame);
                        end
                    end
                end
            else
                for i=1:size(I1,4)
                    for k=1:dm(3)                        
                        subplot(1,3,1);imagesc(rot90(I_ref(:,:,k)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',k)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,3,2);imagesc(rot90(I1(:,:,k,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,3,3);imagesc(rot90(I2(:,:,k,i)),intv_src); axis off; title({'Images 2', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                        pause(Pause)
                        if dummy_movie
                            frame = getframe(gcf);
                            writeVideo(aviobj,frame);
                        end
                    end
                end
            end

        % single slice
        else
            for i=1:size(I1,4)
                subplot(1,3,1);imagesc(rot90(I_ref(:,:,slice)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',slice_for_ref)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                subplot(1,3,2);imagesc(rot90(I1(:,:,slice,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',slice,i)}); pbaspect([x(1) x(2) x(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                subplot(1,3,3);imagesc(rot90(I2(:,:,slice,i)),intv_src); axis off; title({'Images 2', sprintf('Slice %d, Volume %d',slice,i)}); pbaspect([x(1) x(2) x(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                pause(Pause)
                if dummy_movie
                    frame = getframe(gcf);
                    writeVideo(aviobj,frame);
                end
            end
        end

    % 3 datasets
    case 3
        
        % all slices
        if dummy_allslices
            if ~dummy_order
                for k=1:dm(3)
                    for i=1:size(I1,4)
                        subplot(1,4,1);imagesc(rot90(I_ref(:,:,k)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',k)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,4,2);imagesc(rot90(I1(:,:,k,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,4,3);imagesc(rot90(I2(:,:,k,i)),intv_src); axis off; title({'Images 2', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        subplot(1,4,4);imagesc(rot90(I3(:,:,k,i)),intv_src); axis off; title({'Images 3', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k)),N_lines,'r')
                        sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                        pause(Pause)
                        if dummy_movie
                            frame = getframe(gcf);
                            writeVideo(aviobj,frame);
                        end
                    end
                end
            else
                for i=1:size(I1,4)
                    for k=1:dm(3)
                        subplot(1,4,1);imagesc(rot90(I_ref(:,:,k)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',k)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k,i)),N_lines,'r')
                        subplot(1,4,2);imagesc(rot90(I1(:,:,k,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k,i)),N_lines,'r')
                        subplot(1,4,3);imagesc(rot90(I2(:,:,k,i)),intv_src); axis off; title({'Images 2', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k,i)),N_lines,'r')
                        subplot(1,4,4);imagesc(rot90(I3(:,:,k,i)),intv_src); axis off; title({'Images 3', sprintf('Slice %d, Volume %d',k,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                        hold on;
                        contour(rot90(I_contour(:,:,k,i)),N_lines,'r')
                        sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                        pause(Pause)
                        if dummy_movie
                            frame = getframe(gcf);
                            writeVideo(aviobj,frame);
                        end
                    end
                end
            end

        % single slice
        else
            for i=1:size(I1,4)
                subplot(1,4,1);imagesc(rot90(I_ref(:,:,slice)),intv_ref); axis off; title({'Reference image',sprintf('Slice %d',slice)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                subplot(1,4,2);imagesc(rot90(I1(:,:,slice,i)),intv_src); axis off; title({'Images 1', sprintf('Slice %d, Volume %d',slice,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                subplot(1,4,3);imagesc(rot90(I2(:,:,slice,i)),intv_src); axis off; title({'Images 2', sprintf('Slice %d, Volume %d',slice,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                subplot(1,4,4);imagesc(rot90(I3(:,:,slice,i)),intv_src); axis off; title({'Images 3', sprintf('Slice %d, Volume %d',slice,i)}); pbaspect([dm(1) dm(2) dm(3)]);
                hold on;
                contour(rot90(I_contour(:,:,slice)),N_lines,'r')
                sgtitle(contour_text, 'FontSize', 18, 'FontWeight', 'bold')
                pause(Pause)
                if dummy_movie
                    frame = getframe(gcf);
                    writeVideo(aviobj,frame);
                end
            end
        end
end

close(fig1)

if dummy_movie
    close(aviobj);
end
end
function V = acid_load_4Dimage(P)

    if ~isempty(P)
        
        % return error for a single(3D) input
        if size(P,1) == 1
            struct4D = nifti(P);
            dim4D = struct4D.dat.dim;
                 
            if length(dim4D) < 4
                error('A single 3D source image was selected. Choose a single 4D volume or select all 3D volumes manually.');
            else
                n = dim4D(4);
                if n == 1
                    error('A single 3D source image was selected. Choose a single 4D volume or select all 3D volumes manually.');
                end
            end
            
            pos_comma = strfind(P,',');
            if ~isempty(pos_comma)
                P = strcat(repmat(P(1:pos_comma-1), n, 1), ',', num2str([1:n]'));
            else
                P = strcat(repmat(P, n, 1), ',', num2str([1:n]'));
            end
        end
        
        % load in headers
        for i = 1:size(P,1)
            Vtmp = spm_vol(P(i,:));
            if size(Vtmp,1) > 1
                V(i,1) = Vtmp(1);
            else
                V(i,1) = Vtmp;
            end
        end
    else
        V = [];
    end
end
function acid_fusion(P_target,P_b,P_sc,overlap_lim,interp)


% load in target

% P_target=spm_select();
V_target = spm_vol(P_target);



keyword = 'FUSION';

[path,fname,~] = spm_fileparts(V_target{1}.fname);
p_out = acid_bids(path,fname,keyword,1);

% start logging
diary([p_out filesep 'logfile_' keyword '.txt'])


% load in brain image

% P_b=spm_select();
V_b = spm_vol(P_b);
I_b_uncut = acid_read_vols(V_b{1},V_target{1},interp);
I_b = I_b_uncut((overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim);

% load in spinal cord image

% P_sc=spm_select();
V_sc = spm_vol(P_sc);
I_sc_uncut = acid_read_vols(V_sc{1},V_target{1},interp);
I_sc = I_sc_uncut((overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim);


% for i = 1:3
% V_eigenvector_brain = spm_vol(P_eigenvector_brain);
% I_eigenvector_brain_uncut = acid_read_vols(V_eigenvector_brain{i},V_target{1},interp);
% I_eigenvector_brain(:,:,:,i) = I_eigenvector_brain_uncut((overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim);
% end
% 
% for i = 1:3
% V_eigenvector_sc = spm_vol(P_eigenvector_sc);
% I_eigenvector_sc_uncut = acid_read_vols(V_eigenvector_sc{i},V_target{1},interp);
% I_eigenvector_sc(:,:,:,i) = I_eigenvector_sc_uncut((overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim,(overlap_lim+1):end-overlap_lim);
% end


% overlap 
nozeros_Isc = find(I_sc);

indi_nozeros = find(I_b(nozeros_Isc) > 0);

% for i = 1:3
%     I_eigen_out(:,:,:,i) = I_eigenvector_brain(:,:,:,i) + I_eigenvector_sc(:,:,:,i);
% 
%     I_eigen_out(nozeros_Isc(indi_nozeros)) = (I_eigenvector_sc(nozeros_Isc(indi_nozeros)) + I_eigenvector_brain(nozeros_Isc(indi_nozeros))) / 2;
% end

I_out = I_b + I_sc;

I_out(nozeros_Isc(indi_nozeros)) = (I_sc(nozeros_Isc(indi_nozeros)) + I_b(nozeros_Isc(indi_nozeros))) / 2;


acid_write_vol(I_out, V_target{1}, p_out, [keyword '-FA'], 'same', '', '', 1);
% acid_write_vol(I_eigen_out, V_target{1}, p_out, [keyword '-V123'], 'same', '', '', 1);

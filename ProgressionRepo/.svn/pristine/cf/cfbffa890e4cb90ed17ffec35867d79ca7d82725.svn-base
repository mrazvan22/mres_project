function I_new = convert_indimg(I_old, V_old, V_new)

coord_vox_old = zeros(length(I_old), 1);
[coord_vox_old(:, 1), coord_vox_old(:, 2), coord_vox_old(:, 3)] = ...
    ind2sub(V_old.dim, I_old);
coord_vox_old = [coord_vox_old ones(length(I_old), 1)];
coord_vox_new = zeros(size(coord_vox_old));
for vox = 1:length(I_old),
    
    coord_mm = V_old.mat*coord_vox_old(vox, :)';
    coord_vox_new(vox, :) = inv(V_new.mat)*coord_mm;
    
end
coord_vox_new= round(coord_vox_new);
coord_vox_new = unique(coord_vox_new, 'rows');
I_new = sub2ind(V_new.dim, coord_vox_new(:, 1), ...
    coord_vox_new(:, 2), coord_vox_new(:, 3));

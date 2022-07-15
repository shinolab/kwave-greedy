function res = gen_transducers(trans_x, trans_y, trans_spacing)
trans_positions = zeros(trans_x * trans_y, 3);
for iy = 1:trans_y
    for ix = 1:trans_x
        trans_positions(ix + trans_x * (iy-1), :) = [(ix-1) * trans_spacing, (iy-1) * trans_spacing, 0];
    end
end
center = mean(trans_positions);
for iy = 1:trans_y
    for ix = 1:trans_x
        trans_positions(ix + trans_x * (iy-1), :) = trans_positions(ix + trans_x * (iy-1), :) - center;
    end
end
res = trans_positions;
end

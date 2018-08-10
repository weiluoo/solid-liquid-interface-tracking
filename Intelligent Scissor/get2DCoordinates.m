function [row, col] = get2DCoordinates(idx, ncols)

row = mod(idx, ncols) + 1;
col = rem(idx, ncols);
if col == 0
    col = ncols;
end

end
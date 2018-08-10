function p = pathGenerator(startRow, endRow, len)

if startRow ~= endRow
    p = round(startRow:(endRow-startRow)/len:endRow);
else
    p = ones(1, len+1) * startRow;
end

end
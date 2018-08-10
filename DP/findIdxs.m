function idxs = findIdxs(choices, targets)

idxs = [];
for i = 1:size(targets, 1)
    idxs = [idxs; find(choices == targets(i))];
end

end
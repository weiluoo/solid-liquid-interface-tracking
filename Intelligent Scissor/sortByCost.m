function sortedPaths = sortByCost(paths)

[~, idx] = sort(paths(:,1));
sortedPaths = paths(idx,:);

end
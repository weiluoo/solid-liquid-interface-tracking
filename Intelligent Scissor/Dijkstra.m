function shortestPath = Dijkstra(cost, startP, endP)
[rows, cols] = size(cost);
length_penalty = 3.0;


minRow = min(startP(1), endP(1)) - 8;
maxRow = max(startP(1), endP(1)) + 8;
minCol = startP(2);
maxCol = endP(2);
maxLen = sum(abs(endP - startP)) + 8 + 2;   % % % [cost, usage, points' idx] % % %


paths = zeros(1, maxLen);
paths(2) = 3; 
paths(3) = get1DIdx(start(1), start(2), cols);
visited = [];


while true
    
end

end
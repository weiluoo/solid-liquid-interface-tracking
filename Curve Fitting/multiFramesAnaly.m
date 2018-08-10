close all; clear all; clc;
warning off;

contrastPath = '/usr/local/home/wlhz2/Documents/MATLAB/Data/contrastImages/80_Al_CastSubstrate_P100S0.5_U18G15_500ns_S1/';
imageFormat = '*.tif';
imageFiles = dir([contrastPath imageFormat]);

startNum = 178;
endNum = 187;
interval = 3;
n = (endNum - startNum) / interval + 1;
% % % initialize % % % 
image = imread([contrastPath imageFiles(startNum).name]);
% % 
figure(1); imshow(image); hold on;
% % 
Curve = csvread(sprintf('groundTrue/%d.csv', startNum));
y = Curve(:,6); x = Curve(:,7);
Poly = createFit(y, x);
Vals = coeffvalues(Poly);
Derivative = Vals(1:9) .* [9:-1:1];
Y = [min(y):0.1:max(y)]';
X = Poly(Y);
% % 
plot(Y, X, 'Color', 'w', 'LineWidth', 1)
% % 

N = size(Y, 1);
Ps = zeros(N, 2*n);
Directs = zeros(N, n);
Dists = zeros(N, n-1);
Ps(:,1:2) = [X, Y];
count = 0;

for jNum = (startNum+interval):interval:endNum
    Curve = csvread(sprintf('groundTrue/%d.csv', jNum));
    y = Curve(:,6); x = Curve(:,7);
    Poly = createFit(y, x);
    jY = [min(y):0.1:max(y)]';
    jX = Poly(jY);
    jN = length(jY);
    
    % % % Normal lines for previous solid-liquid interface % % %
    idxs = (Ps(:, count*2+1) ~= 0) & (Ps(:, count*2+2) ~= 0);
    iN = sum(idxs);
    iX = Ps(idxs, count*2+1); iY = Ps(idxs, count*2+2);
    ias = -1 ./ Poly8(Derivative, iY);
    ibs = iX - ias .* iY;
    NormalVals = repmat(ias, [1, jN]) .* repmat(jY', [iN, 1]) + repmat(ibs, [1, jN]);
    
    % % 
    plot(jY, jX, 'Color', 'y', 'LineWidth', 1)
    % % 
    Vals = coeffvalues(Poly);
    Derivative = Vals(1:9) .* [9:-1:1];
    
    % % % Find corresponding points in current interface % % %
    crossDists = repmat(jX', [iN, 1]) - NormalVals;
    overZeros = crossDists >= 0;
    valids = ones(iN, 1);
    tmp = sum(overZeros, 2);
    valids = valids - (tmp == 0 | tmp == jN);
    [~, minIdxs] = min(abs(crossDists), [], 2);
    jX = jX(minIdxs); jY = jY(minIdxs);
    
    % % % %
    for i = 1:50:iN
        if valids(i) == 0
            continue;
        end
        plot([jY(i); iY(i)], [jX(i); iX(i)], 'Color', 'b', 'LineWidth', 0.1);
    end
    % % % %
    
    % % % Compute the speed % % %
    dists = sqrt(((iX - jX) .^ 2) + ((iY - jY) .^ 2));

    % % %
    count = count + 1;
    Ps(idxs, count*2+1) = jX .* valids;
    Ps(idxs, count*2+2) = jY .* valids;
    Dists(idxs, count) = dists .* valids;
    Directs(idxs, count) = ias;

end

figure(2); imshow(image); hold on;
for i = 1:50:N
    x = Ps(i, 1); y = Ps(i, 2);
    for j = 1:n-1
        if (Ps(i, j*2+1) ~= 0) & (Ps(i, j*2+2) ~= 0)
            x = [x; Ps(i, j*2+1)]; y = [y; Ps(i, j*2+2)];
        else
            break;
        end
    end
    plot(y, x, 'Color', 'g', 'LineWidth', 0.1);
end
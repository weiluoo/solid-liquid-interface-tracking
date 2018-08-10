close all; clear all; clc;

load('allTraces.mat', 'allTraces');

iNum = 70;
jNum = 72;

iTrace = allTraces(:,:,iNum);
jTrace = allTraces(:,:,jNum);

[ix, iy] = find(iTrace);
[jx, jy] = find(jTrace);

% iN = length(iy);
% jN = length(jy);
% 
% iPoint = [ix(1), iy(1)];
% iPoints = iPoint;
% for i = 2:iN
%     if (ix(i) ~= iPoint(1))
%         iPoint = [ix(i), iy(i)];
%         iPoints = [iPoints; iPoint];
%     end
% end
% 
% jPoint = [jx(1), jy(1)];
% jPoints = jPoint;
% for j = 2:jN
%     if (jx(j) ~= jPoint(1))
%         jPoint = [jx(j), jy(j)];
%         jPoints = [jPoints; jPoint];
%     end
% end
% 
% figure(1); imshow(iTrace); title(num2str(iNum)); hold on;
% plot(iPoints(:,2), iPoints(:,1), 'yx');
% figure(2); imshow(jTrace); title(num2str(jNum)); hold on;
% plot(jPoints(:,2), jPoints(:,1), 'yx');
% 
% iPoly = createFit(iPoints(:,2), iPoints(:,1));
% jPoly = createFit(jPoints(:,2), jPoints(:,1));

iPoly = createFit(iy, ix);
jPoly = createFit(jy, jx);

iVals = coeffvalues(iPoly);
iDerivative = iVals(1:9) .* [9:-1:1];

% % % visualize % % % 
% iPoints = [iPoly(iy(1:30:end)), iy(1:30:end)];
% jPoints = [jPoly(jy(1:30:end)), jy(1:30:end)];
figure(1); imshow(ones(size(iTrace)), []); hold on;
plot(iPoly, 'r-', [iy(1); iy(end)], [iPoly(iy(1)); iPoly(iy(end))], 'k.');
plot(jPoly, 'b-', [jy(1); jy(end)], [jPoly(jy(1)); jPoly(jy(end))], 'k.');
% % % % % % % % % % % 

% % % common column % % % 
commonY = [max(min(iy), min(jy)):0.1:min(max(iy), max(jy))]';
iXs = iPoly(commonY);
jXs = jPoly(commonY);
N = length(commonY);

figure(1); imshow(image); hold on;
plot(commonY, iXs, 'Color', 'b', 'LineWidth', 1)
plot(commonY, iNormalVals(3611,:), 'Color', 'k', 'LineWidth', 1)
% % % % % % % % % % % % % 

% % % compute the normal of line i % % % 
iValsPoly8 = Poly8(iDerivative, commonY);
ias = -1 ./ iValsPoly8;
ibs = iXs - ias .* commonY;
iNormalVals = repmat(ias, [1, N]) .* repmat(commonY', [N, 1]) + repmat(ibs, [1, N]);
crossDists = repmat(jXs', [N, 1]) - iNormalVals;
overZeros = crossDists >= 0;
% % % % % % % % % % % % % % % % % % % % % 

% % % find the cross points % % % 
valids = ones(N, 1);
tmp = sum(overZeros, 2);
valids = valids - (tmp == 0 | tmp == N);
[~, idx] = min(abs(crossDists), [], 2);
crossPoints = [jXs(idx), commonY(idx)];
% figure(1); hold on;
% plot(crossPoints(:,2), crossPoints(:,1), 'gx');
% % % % % % % % % % % % % % % % % 

% % % compute the distance % % %
dists = sqrt(((crossPoints(:,1) - iXs) .^ 2) + ((crossPoints(:,2) - commonY) .^ 2));
dists = dists .* valids;
directs = (crossPoints(:,1) - iXs) > 0;
% % % % % % % % % % % % % % % % % 

figure(1); imshow(ones(size(iTrace)), []); hold on;
plot(commonY, iXs, 'Color', 'b', 'LineWidth', 0.1);
for i = 1:10:N
    if valids(i) == 0
        continue;
    end
    
    if directs(i) == 1
        plot([commonY(i); crossPoints(i, 2)], [iXs(i); crossPoints(i, 1)], 'Color', 'r', 'LineWidth', 0.1);
    else
        plot([commonY(i); crossPoints(i, 2)], [iXs(i); crossPoints(i, 1)], 'Color', 'g', 'LineWidth', 0.1)
    end
end




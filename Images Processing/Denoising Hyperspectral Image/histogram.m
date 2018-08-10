clc; %clear the command window 
close all; %close all figure windows
clear all; %clear all variables in the workspace

A = imread('1.jpg');
B = imread('2.jpg');
C = imread('3.jpg');
D = imread('5.jpg');
E = imread('6.jpg');
F = imread('8.jpg');
G = imread('9.jpg');
imhist(A);title('original image');
figure;
imhist(B);title('restored image with method 1');
figure;
imhist(C);title('noise image with method 1');
figure;
imhist(D);title('restored image with method 2');
figure;
imhist(E);title('noise image with method 2');
figure;
imhist(F);title('restored image with method 3');
figure;
imhist(G);title('noise image with method 3');


% A = imread('57.jpg');
% B = imread('65.jpg');
% imwrite([A',B'],[num2str(73) '.jpg']);
% 
% C = imread('60.jpg');
% D = imread('68.jpg');
% imwrite([C',D'],[num2str(74) '.jpg']);
% 
% G = imread('61.jpg');
% H = imread('69.jpg');
% imwrite([G',H'],[num2str(75) '.jpg']);
% 
% E = imread('62.jpg');
% F = imread('70.jpg');
% imwrite([E',F'],[num2str(76) '.jpg']);
% 
% I = imread('63.jpg');
% J = imread('71.jpg');
% imwrite([I',J'],[num2str(77) '.jpg']);
% 
% K = imread('64.jpg');
% L = imread('72.jpg');
% imwrite([K',L'],[num2str(78) '.jpg']);
% imread('dc.tif');
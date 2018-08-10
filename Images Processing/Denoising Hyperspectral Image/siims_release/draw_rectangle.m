function uw = draw_rectangle( u,x1,x2,y1,y2,b,c)
%DRAW_RECTANGLE 
%
%   uw = draw_rectangle( u,x1,x2,y1,y2,b,c)
%
%   Draw a rectangle on the image u (the content of the image on the
%   rectangle border is overwritten).
%
%   INPUT
%   u           Input image
%   x1,y1       Up left coordinates of the rectangle
%   x2,y2       Down right coordinates of the rectangle
%   b           Width of the rectangle border
%   c           Color of the rectangle
%
%   OUTPUT
%   uw          Output image

if nargin<7
    c = [1 0 0];
end
c = reshape(c,[1 1 3]);

if size(u,3)>1
    uw = u;
else
    uw = repmat(u,[1 1 3]);
end
uw(y1:y1+b-1,x1:x2,:) = repmat(c,[b,x2-x1+1,1]);
uw(y2-b+1:y2,x1:x2,:) = repmat(c,[b,x2-x1+1,1]);
uw(y1:y2,x1:x1+b-1,:) = repmat(c,[y2-y1+1,b,1]);
uw(y1:y2,x2-b+1:x2,:) = repmat(c,[y2-y1+1,b,1]);


end


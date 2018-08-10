function d = disc( M,N,cx,cy,r,v )
%DISC Compute a disc image
%   d = disc( M,N,cx,cy,r,v ) with parameters 
%   [M,N] : size of the output image
%   [cx,cy] : center of the disc
%   r : disc radius
%   v : maximum value of the image

[Y,X] = meshgrid(1:N,1:M);
d = zeros(M,N);
d((X-cx).^2+(Y-cy).^2<=(r+0.5)^2) = v;
end


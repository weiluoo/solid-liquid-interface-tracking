function co = get_conditioning_points( mask,w )
%GET_CONDITIONING_POINTS FROM THE MASK
%   c = get_conditioning_points( mask,w )
%   Image c is 1 on a border of the mask with width w.
%   mask is a binary image (indicator function of a mask)
%
%   NB: the mask is supposed to be the same on all color channels
%
%   Author: Arthur Leclaire
%   v 1.0 (03/2016) First public version
    
    C = size(mask,3);
    if C>1
        mask = mask(:,:,1);
    end
    [M,N] = size(mask);
    M2 = 2*M; N2 = 2*N;
    bmask = zeros(M2,N2);
    bmask(1:M,1:N) = double(mask>1e-10);
    k = zeros(M2,N2);
    k(1:2*w+1,1:2*w+1) = 1;
    k = circshift(k,[-w -w]);
    co = real(ifft2(fft2(bmask).*fft2(k)));
    co = (co>0.5);
    co(bmask>0) = 0;
    co = co(1:M,1:N);
    co = repmat(co,[1 1 C]);
end


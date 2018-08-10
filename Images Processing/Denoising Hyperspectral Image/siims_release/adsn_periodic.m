function out = adsn_periodic( s , mu )
%ADSN_PERIODIC Compute a periodic Asymptotic Discrete Spot Noise texture.
%
%   out = adsn_periodic(s,mu) 
%
%   computes a realization of the Gaussian 
%   circularly stationary random field of mean mu and whose 
%   covariance function is the periodic autocorrelation of s.
%
%   Notice that the covariance of the resulting field is the
%   periodic autocorrelation of s.
%
%   NB : 
%   - The mean value of s is not substracted.
%   - If size(s,1)>M or size(s,2)>N , then s is cropped.
%   - The input s can have multiple channels.
%
%   This texture model is presented in the paper
%       "Random Phase Textures: Theory and Synthesis", 
%       (B. Galerne, Y. Gousseau, J.-M. Morel), 
%       IEEE Transactions on Image Processing, 2011.
%
%   Author : Arthur Leclaire
%   v. 1.0 (03/2014) : first public version

if (size(s,3)~=size(mu,3))
    error('s and mu must have the same number of channels')
end

[M,N,C] = size(mu);

if size(s,1)>M
    s = s(1:M,:,:);
end
if size(s,2)>N
    s = s(:,1:N,:);
end
m = size(s,1);
n = size(s,2);

% Convolution of the kernel s with a normal Gaussian white noise.
out = zeros(M,N,C);
W = randn(M,N);
fW = fft2(W);
for c=1:C
    out(:,:,c) = [s(:,:,c) zeros(m,N-n); zeros(M-m,N)];
    out(:,:,c) = real(ifft2(fft2(out(:,:,c)).*fW));
end
% Reimpose the mean image
out = mu + out;

end
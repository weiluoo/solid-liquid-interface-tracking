function v = gaussian_convol( u,s )
%GAUSSIAN_CONVOL Gaussian periodic convolution

    [M,N] = size(u);
    [Y,X] = meshgrid(0:N-1,0:M-1);
    X(X>M/2) = X(X>M/2)-M;
    Y(Y>N/2) = Y(Y>N/2)-N;
    fk = exp(-2*(s*pi*X/M).^2 - 2*(s*pi*Y/N).^2);
    v = real(ifft2(fk.*fft2(u)));
    

end


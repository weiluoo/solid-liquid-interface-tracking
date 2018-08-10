function v = mtimes(A,u)
% Apply the covariance operator A to an image u.

[M,N,C] = size(u);
tmp = zeros(M,N);
v = zeros(M,N,C);
ft = fft2(A.t);
for c=1:C
    tmp = tmp + conj(ft(:,:,c)).*fft2(u(:,:,c));
end
for c=1:C
    v(:,:,c) = real(ifft2(ft(:,:,c).*tmp));
end
v = v + (A.reg)^2 * u;

v(A.ind==0) = 0;

function v = zeropad( u,M,N )
%ZEROPAD v = zeropad( u,M,N )
%   Extend u by zero on a domain of size M x N

[m,n,C] = size(u);
if (m>M)||(n>N)
    warning('Extended domain smaller than domain of u')
end

v = zeros(M,N,C);
v(1:min(M,m),1:min(N,n),:) = u(1:min(M,m),1:min(N,n),:);

end


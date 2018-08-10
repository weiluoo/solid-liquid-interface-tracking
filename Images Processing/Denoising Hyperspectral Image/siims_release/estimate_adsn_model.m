function [ t,m ] = estimate_adsn_model( u,Mb,Nb,maskrgb )
%ESTIMATE_GAUSSIAN_MODEL COMPUTE THE MEAN AND TEXTON ASSOCIATED TO U
%   [ t,m ] = estimate_adsn_model( u,Mb,Nb,maskrgb )
%   
%   INPUT
%   u       Original texture image
%   Mb,Nb   (optional) The texton is embedded in a MbxNb image
%   maskrgb (optional) Estimate the ADSN model outside a mask
%       
%   OUTPUT
%   t       Texton of the ADSN model
%   m       Mean value of the ADSN model
%
%   Author: Arthur Leclaire
%   v 1.0 (03/2016) First public version

[M,N,C] = size(u);
if nargin<2
    Mb = M; Nb = N;
end

if nargin<4
    m = mean(mean(u,2));
    tl = 1/sqrt(M*N)*(u-repmat(m,[M N]));
    
    t = zeros(Mb,Nb,C);
    t(1:M,1:N,:) = tl(1:M,1:N,:);
else
    m = zeros(1,1,C);
    t = zeros(Mb,Nb,C);
    maskrgb = (abs(maskrgb)>1e-10);
    k = sqrt(C/nnz(~maskrgb));
    for c=1:C
        uc = u(:,:,c); mc = maskrgb(:,:,c);
        m(1,1,c) = mean(mean(uc(mc==0)));
        tc = zeros(M,N,1);
        tc(mc==0) = k*(uc(mc==0)-m(1,1,c));
        t(1:M,1:N,c) = tc;
    end
end

end


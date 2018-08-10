function [kcomp,condcov,kcoeff,kcompvisu,v,GammaMix,GammaCond,GammaMask,GammaCondPinv] = get_conditional_model(u,mv,t,indm,indc,options)
%   [kcomp,condcov,kcoeff,kcompvisu,v,GammaMix,GammaCond,GammaMask,GammaCondPinv] = get_conditional_model(u,mv,t,indm,indc,options)
%
%   INPUT
%   u       Initial image
%   mv      Mean value of the estimated Gaussian model (ADSN)
%   t       Texton of the estimated Gaussian model (ADSN)
%   indm    Indicator function of masked region
%   indc    Indicator function of conditioning region
%   options Structure containing optional inputs:
%   - options.per  (default 0)
%       Use periodic ADSN if per==1, otherwise use non-periodic ADSN
%   - options.reg (default 0)
%       Add regularization in the linear system
%
%   OUTPUT

%   NB:
%   - The mean value is reimposed on kvisu and innov for visualization.
%   - indm can have dimensions MxNxC (provided that the
%      arrays indm(:,:,c) are the same for c=1:C). Idem for indc.
%

options.trash = 0;

if (nargin<6)||(~isfield(options,'per'))
    per = 0;
else
    per = options.per;
end

if (nargin<6)||(~isfield(options,'reg'))
    reg = 0;
else
    reg = options.reg;
end

[M,N,C] = size(u);

% Handle the dimensions of indm and indc
indm = indm(:,:,1);
indc = indc(:,:,1);

% technical restriction
% (in order to simplify the covariance computation)
if (size(t,1)>M)||(size(t,2)>N)
    error('The texton must not be larger than the input');
end

if per==0
    M = 2*M; N = 2*N;
    u = zeropad(u,M,N);
    t = zeropad(t,M,N);
    indm = zeropad(indm,M,N);
    indc = zeropad(indc,M,N);
end

mp = find(indm>0);
mn = length(mp);
[ma,mb] = ind2sub([M N],mp);
cp = find(indc>0);
cn = length(cp);
[ca,cb] = ind2sub([M N],cp);

% Covariance the whole random field
covt = zeros(M,N,C,C);
ft = fft2(t);
for a=1:C
    for b=1:C
        covt(:,:,a,b) = real(ifft2(conj(ft(:,:,a)).*ft(:,:,b)));
    end
end

% Fill the matrix of the system
GammaCond = zeros(C*cn);
for j=1:cn
    for a=1:C
        for b=1:C
            shifta = cn*(a-1);
            shiftb = cn*(b-1);
            % covariance around (ca(j),cb(j))
            covtt = circshift(covt(:,:,a,b),[ca(j)-1 cb(j)-1]);
            % extract the part on the conditioning points
            GammaCond(shifta+j,shiftb+(1:cn)) = covtt(cp);
        end
    end
end
GammaCond = GammaCond + reg^2*eye(C*cn);

% Fill the matrix of the system
GammaMask = zeros(C*mn);
for j=1:mn
    for a=1:C
        for b=1:C
            shifta = mn*(a-1);
            shiftb = mn*(b-1);
            % covariance around (ma(j),mb(j))
            covtt = circshift(covt(:,:,a,b),[ma(j)-1 mb(j)-1]);
            % extract the part on the conditioning points
            GammaMask(shifta+j,shiftb+(1:mn)) = covtt(mp);
        end
    end
end
GammaCond = GammaCond + reg^2*eye(C*cn);

GammaMix = zeros(C*cn,C*mn);
for j=1:cn
    for a=1:C
        for b=1:C
            shifta = cn*(a-1);
            shiftb = mn*(b-1);
            % covariance around (ca(j),cb(j))
            covtt = circshift(covt(:,:,a,b),[ca(j)-1 cb(j)-1]);
            % extract the part on the conditioning points
            GammaMix(shifta+j,shiftb+(1:mn)) = covtt(mp);
        end
    end
end

disp(['pinv tolerance value = ' num2str(max(size(GammaCond))*norm(GammaCond)*eps)]);
GammaCondPinv = pinv(GammaCond);

kcoeff = GammaMix'*GammaCondPinv;

% draw an independent realization of the random field
z = adsn_periodic(t,zeros(M,N,C))+reg*randn(M,N,C);

% Compute centered conditioning values
condval = zeros(C*cn,1);
zval = zeros(C*cn,1);
for a=1:C
    shifta = cn*(a-1);
    uc = u(:,:,a);
    zc = z(:,:,a);
    condval(shifta+(1:cn)) = uc(cp)-mv(:,:,a);
    zval(shifta+(1:cn)) = zc(cp);
end

% Compute kriging component
kcomp = kcoeff*condval;
kcompnoise = kcoeff*zval;

kcompvisu = u;
v = u;
for a=1:C
    shifta = mn*(a-1);
    uc = u(:,:,a);
    uc(mp) = mv(:,:,a)+kcomp(shifta+(1:mn));
    kcompvisu(:,:,a) = uc;
    vc = v(:,:,a); zc = z(:,:,a);
    vc(mp) = uc(mp) + zc(mp) - kcompnoise(shifta+(1:mn));
    v(:,:,a) = vc;
end

%condcov = GammaMask - kcoeff*GammaMix;

 condcov = GammaMask - kcoeff*GammaMix - GammaMix'*kcoeff' + ...
   kcoeff*GammaCond*kcoeff';

if per==0
    M = floor(M/2); N = floor(N/2);
    kcompvisu = kcompvisu(1:M,1:N,:);
    v = v(1:M,1:N,:);
end
end

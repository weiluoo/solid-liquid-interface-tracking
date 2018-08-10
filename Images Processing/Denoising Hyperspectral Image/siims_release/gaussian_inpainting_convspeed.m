function [v,distvref,distiter,resl2norm,reslinfnorm] = gaussian_inpainting_convspeed(u,mv,t,indm,indc,options)
%GAUSSIAN_INPAINTING_CONVSPEED
% [v,kvisu,innov,] = gaussian_inpainting_convspeed(u,mv,t,indm,indc,options)
%
%   THIS FUNCTION WAS DERIVED FROM gaussian_inpainting (version of the
%   19 april 2016) IN ORDER TO ANALYZE THE CONVERGENCE SPEED OF THE ALGO.
%
%   INPUT
%   u       Initial image
%   mv      Mean value of the estimated Gaussian model (ADSN)
%   t       Texton of the estimated Gaussian model (ADSN)
%   indm    Indicator function of masked region (same size as u)
%   indc    Indicator function of conditioning region (same size as u)
%   options Structure containing optional inputs:
%   - options.per  (default 0) 
%       Use periodic ADSN if per==1, otherwise use non-periodic ADSN
%   - options.ep (default 1e-3)
%       Desired precision for the solution of the linear system solved
%       by the conjugate gradient descent
%   - options.imax (default 1000)
%       Maximum number of iterations in gradient descent
%   - options.normal (default 1)
%       Perform conjugate gradient descent on the normal equations.
%   - options.verb (default 0)
%       Verbose mode
%   - options.vref
%       Reference image to analyze the convergence of the algo.
%   - options.reg (default 0)
%       Add regularization in the linear system
%
%   OUTPUT
%   v         Inpainted image
%   distvref    Distance to the reference image (along with the iterations)
%   distiter    Distance between two iterates
%   resl2norm   Residual L2-norm
%      (normalized by sqrt(number of conditioning points x nb channels))
%   reslinfnorm Residual Linf-norm
%
%   NB: 
%   - The mean value is reimposed on kvisu and innov for visualization.
%   - indm can have dimensions MxNxC (provided that the 
%      arrays indm(:,:,c) are the same for c=1:C). Idem for indc.
%

options.trash = 0;

if ~isfield(options,'per')
    per = 0;
else
    per = options.per;
end

if ~isfield(options,'ep')
    ep = 1e-3;
else
    ep = options.ep;
end

if ~isfield(options,'imax')
    imax = 1000; % max number of iterations
else
    imax = options.imax;
end

if ~isfield(options,'reg')
    reg = 0;
else
    reg = options.reg;
end

% In this code, we always use the normal equations.
%
% if ~isfield(options,'normal')
%     normal = 1;
% else
%     normal = options.normal;
% end

if ~isfield(options,'verb')
    verb = 0;
else
    verb = options.verb;
end

[M,N,C] = size(u);

% Handle the dimensions of indm and indc
indm = repmat(indm(:,:,1),[1 1 C]);
indc = repmat(indc(:,:,1),[1 1 C]);    
nc = sum(indc(:)); % nb cond points x nb channels

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

indm = logical(indm);
indc = logical(indc);

mv = repmat(mv,[M N 1]);
v = u;

G = convcovr(t,ones(M,N,C),reg); % convolution by the covariance
A = convcovr(t,indc,reg); % convolution by the cov and restriction to indc

% draw an independent realization of the random field
z = adsn_periodic(t,zeros(M,N,C)) + reg*randn(M,N,C);
% NB: if per==1, t has already been extended by zero-padding
% so that the half-crop of z is indeed a sample of the
% non-periodic ADSN associated to t.

% if kvisu and innov are not required (i.e. nargout==1),
% then one only needs the kriging component of u-mv-z
%   (faster to compute)

rhs = zeros(M,N,C); rhs(indc) = u(indc)-mv(indc)-z(indc);

distvref = zeros(imax,1);
distiter = zeros(imax,1);
resl2norm = zeros(imax,1);
reslinfnorm = zeros(imax,1);

% initialization
i = 0;
x = zeros(M,N,C);
r = A'*rhs; 
nr2 = sum(r(:).^2);

p = r;
rninf = max(abs(r(:))); % L^\infty norm of the residual
stop = 0;
% main loop
while i<imax && rninf>ep && ~stop
    i = i+1;
    AtAp = A'*(A*p);
    spAtAp = sum(p(:).*AtAp(:));
    alpha = nr2/spAtAp;
    x = x + alpha*p;
    r = r - alpha*AtAp;
    nrold2 = nr2;
    nr2 = sum(r(:).^2);
    rninf = max(abs(r(:)));
    
    beta = nr2/nrold2;
    p = r + beta*p;
    
    resl2norm(i) = sqrt(nr2/nc);
    reslinfnorm(i) = rninf;
    
    if verb==1
        disp(['iteration ' num2str(i) ...
            ', Residual inf norm = ' num2str(rninf) ...
            ', Residual l2 norm = ' num2str(sqrt(nr2/nc))])
    elseif verb==2
        % print also the value of the functional
        value = sum(sum(sum((A*x-rhs).^2)));
        disp(['iteration ' num2str(i) ...
            ', Residual inf norm = ' num2str(rninf) ...
            ', Value = ' num2str(value)]);
    end
    
    kp = G*x;
    vold = v;
    v(indm) = mv(indm) + kp(indm) + z(indm);

    if isfield(options,'visu')
        imshow(v(1:floor(M/2),1:floor(N/2),:))
        drawnow
    end
    distiter(i) = norm(v(:)-vold(:));
    if isfield(options,'vref')
        if per
            vp = v;
        else
            vp = v(1:floor(M/2),1:floor(N/2),:);
        end
        distvref(i) = norm(vp(:)-options.vref(:));
    end
end    

if per==0
    M = floor(M/2); N = floor(N/2);
    v = v(1:M,1:N,:);
end
end

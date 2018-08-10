function [v,kvisu,innov,precsys] = gaussian_inpainting(u,mv,t,indm,indc,options)
%GAUSSIAN_INPAINTING
% [v,kvisu,innov,precsys] = gaussian_inpainting(u,mv,t,indm,indc,options)
%  Microtexture Inpainting by Conditional Gaussian Simulation
%  This algorithm consists in
%   - estimation of a Gaussian model on the masked exemplar,
%   - conditional simulation of this model knowing the values
%       on a set of conditioning points.
%  The linear system involved in conditional simulation is solved
%   using a conjugate gradient descent algorithm.
%   
%   See detailed description of the algorithm in the paper
%   "Gaussian Texture Inpainting" 
%   (Bruno Galerne, Arthur Leclaire)
%   Preprint MAP5 n°2016-25
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
%   - options.ep (default 1e-3 )
%       Desired precision for the solution of the linear system solved
%       by the conjugate gradient descent
%   - options.imax (default 1000)
%       Maximum number of iterations in gradient descent
%   - options.normal (default 1)
%       Perform conjugate gradient descent on the normal equations.
%   - options.verb (default 0)
%       Verbose mode
%   - options.reg (default 0)
%       Add regularization in the linear system
%
%   OUTPUT
%   v       Inpainted image
%   kvisu   (optional) Kriging component
%   innov   (optional) Innovation component
%   precsys (optional) Precision on the solution
%
%   NB: 
%   - The mean value is reimposed on kvisu and innov for visualization.
%   - indm can have dimensions MxNxC (provided that the 
%      arrays indm(:,:,c) are the same for c=1:C). Idem for indc.
%
%   Author: Arthur Leclaire
%   v 1.0 (12/2016) First public version

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

if ~isfield(options,'normal')
    normal = 1;
else
    normal = options.normal;
end

if normal
    disp(['Conjugate gradient descent on normal equations'])
else
    disp('Conjugate gradient descent')
end

if ~isfield(options,'verb')
    verb = 0;
else
    verb = options.verb;
end

[M,N,C] = size(u);

% Handle the dimensions of indm and indc
indm = repmat(indm(:,:,1),[1 1 C]);
indc = repmat(indc(:,:,1),[1 1 C]);    

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
if nargout>1
    
    rhs = zeros(M,N,C); rhs(indc) = u(indc)-mv(indc);
    if normal
        [k,precsys] = solvecgdnormal(A,rhs,zeros(M,N,C),ep,imax,verb);
    else
        [k,precsys] = solvecgd(A,rhs,zeros(M,N,C),ep,imax,verb);
    end
    
    rhs = zeros(M,N,C); rhs(indc) = z(indc);
    if normal
        zk = solvecgdnormal(A,rhs,zeros(M,N,C),ep,imax,verb);
    else
        zk = solvecgd(A,rhs,zeros(M,N,C),ep,imax,verb);
    end
    
    % convolution by the covariance
    k = G*k; 
    zk = G*zk;
    
    v(indm) = mv(indm) + k(indm) + z(indm) - zk(indm);
    kvisu = u;
    kvisu(indm) = mv(indm) + k(indm);
    innov = mv;
    innov(indm) = mv(indm) + z(indm) - zk(indm);

else
    
    rhs = zeros(M,N,C); rhs(indc) = u(indc)-mv(indc)-z(indc);
    if normal
        [k,precsys] = solvecgdnormal(A,rhs,zeros(M,N,C),ep,imax,verb);
    else
        [k,precsys] = solvecgd(A,rhs,zeros(M,N,C),ep,imax,verb);
    end
    k = G*k;
    
    v(indm) = mv(indm) + k(indm) + z(indm);
    
end

if per==0
    M = floor(M/2); N = floor(N/2);
    v = v(1:M,1:N,:);
    if nargout>1
        kvisu = kvisu(1:M,1:N,:);
        innov = innov(1:M,1:N,:);
    end
end
end

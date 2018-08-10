function  res = convcovr(t,ind,reg)
% Define the covariance operator associated to the texton t.
% The operator can be further resctricted to a domain (putting zero
% outside), and also regularized.
%
% INPUT
% t 	texton associated to the Gaussian model
% ind 	indicator function of the restriction region
% reg 	regularization parameter
%
% OUTPUT
% res   covariance operator

res.t = t;
if nargin < 2
  res.ind = ones(size(t));
else
  res.ind = ind;
end
if nargin < 3
    res.reg = 0;
else
    res.reg = reg;
end

res = class(res,'convcovr');

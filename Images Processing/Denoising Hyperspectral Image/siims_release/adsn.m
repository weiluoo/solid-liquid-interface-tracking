function out = adsn( s , mu )
%ADSN Compute an Asymptotic Discrete Spot Noise texture.
%
%   out = adsn(s,mu) computes a realization of the Gaussian stationary 
%   random field of mean mu and whose covariance function is 
%   the autocorrelation of s.
%   The output is of same size as the mean image mu.
%
%   Notice that the covariance of the resulting field is 
%   the non-periodic autocorrelation of s.
%
%   NB : 
%   - The mean value of s is not substracted.
%   - The input s can have multiple channels.
%
%   Author : Arthur Leclaire
%   v. 1.0 (03/2014) : first public version

[M,N,C] = size(mu);
[m,n,~] = size(s);

out = adsn_periodic(s,zeros(M+m,N+n,C));
out = mu + out(1:M,1:N,:);

end
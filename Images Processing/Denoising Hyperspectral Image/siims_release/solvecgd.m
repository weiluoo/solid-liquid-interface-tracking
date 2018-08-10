
function [x,rninf] = solvecgd(A,y,init,ep,imax,verb)
% [x,rninf] = solvecgd(A,y,init,ep,imax,verb)
% Solve Ax = y (A sym def pos) with conjugate gradient descent
%
% Input :
%   A : system matrix
%   y : right-hand side
%   init : initialization of the descent algorithm
%   ep (optional) : desired precision on the L2 norm of the residual Ax-y
%   imax (optional) : maximum number of iterations
%   verb (optional) : 
%    if verb==1, print the residual norm
%    if verb==2, print also the values of
%           the function 1/2 <Ax,x> - <y,x> 
% Output
%   x : system solution
%   rninf : L_\infty norm of the residual Ax-y

if nargin<5
    imax = Inf;
end
if nargin<6
    verb = 0;
end

tic;

i = 0;
% initialization
x = init;
r = y;
nr2 = sum(r(:).^2);

p = r;
rninf = max(abs(r(:))); % L^\infty norm of the residual
rn2 = sqrt(nr2);
stop = 0;
% main loop
while i<imax && rn2>ep && ~stop
    i = i+1;
    Ap = A*p;
    spAp = sum(p(:).*Ap(:));
    alpha = nr2/spAp;
    x = x + alpha*p;
    r = r - alpha*Ap;
    nrold2 = nr2;
    nr2 = sum(r(:).^2);
    rninf = max(abs(r(:)));
    rn2 = sqrt(nr2);
    
    beta = nr2/nrold2;
    p = r + beta*p;
    
    if verb==1
        disp(['iteration ' num2str(i) ...
            ', Residual L2 norm = ' num2str(rn2)]);
    elseif verb==2
        % print also the value of the functional
        value = sum(sum(sum((1/2*(A*x)-y).*x)));
        disp(['iteration ' num2str(i) ...
            ', Residual L2 norm = ' num2str(rn2) ...
            ', Value = ' num2str(value)]);
    end
end

disp(['Time to solve the linear system : ' num2str(toc)]);

end

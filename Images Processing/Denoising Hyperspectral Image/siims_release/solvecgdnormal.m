
function [x,rninf] = solvecgdnormal(A,y,init,ep,imax,verb)
% [x,rninf] = solvecgdnormal(A,y,init,ep,imax,verb)
% Solve Ax = y with conjugate gradient descent
%   on the normal equation.
% This is equivalent to the conjugate gradient descent algorithm to
%   minimize the function |Ax-y|^2 .
%
% Input :
%   A : system matrix
%   y : right-hand side
%   init : initialization of the descent algorithm
%   ep (optional) : desired precision on the L2 norm of the residual Ax-y
%   imax (optional) : maximum number of iterations
%   verb (optional) : 
%       if verb==1, print the residual norm
%   	if verb==2, print also the values of
%           the function |Ax-y|^2
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
r = A'*y; 
nr2 = sum(r(:).^2);

p = r;
rninf = max(abs(r(:))); % L^\infty norm of the residual
rn2 = sqrt(nr2);
stop = 0;
% main loop
while i<imax && rn2>ep && ~stop
    i = i+1;
    AtAp = A'*(A*p);
    spAtAp = sum(p(:).*AtAp(:));
    alpha = nr2/spAtAp;
    x = x + alpha*p;
    r = r - alpha*AtAp;
    nrold2 = nr2;
    nr2 = sum(r(:).^2);
    rninf = max(abs(r(:)));
    rn2 = sqrt(nr2);
    
    beta = nr2/nrold2;
    p = r + beta*p;
    
    if verb==1
        disp(['iteration ' num2str(i) ...
            ', L2 Residual norm = ' num2str(rn2)])
    elseif verb==2
        % print also the value of the functional
        value = sum(sum(sum((A*x-y).^2)));
        disp(['iteration ' num2str(i) ...
            ', L2 Residual norm = ' num2str(rn2) ...
            ', Value = ' num2str(value)]);
    end
end

disp(['Time to solve the linear system : ' num2str(toc)]);

end

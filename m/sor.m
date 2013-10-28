% ------------------------------------------------------------------------------
% FUNCTION:
%       sor
%
% PARAMS:
%       A - <nxn> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el sistema Ax=b mediante el metodo SOR. No siempre hay 
%       convergencia.
% ------------------------------------------------------------------------------

function x = sor(A,b)
    [m n] = size(A);
    tol  = 0.0001;
    w = 1.25;
    xold = ones(n,1)*10000;
    xnew = ones(n,1);
    while ( norm(xnew-xold,"inf")/norm(xnew,"inf") > tol )
        xold = xnew;
        for i = 1:n
            xnew(i) = (w/A(i,i))*(
                        b(i)
                      - sum(A(i,1:i-1)*xnew(1:i-1))
                      - sum(A(i,i+1:n)*xold(i+1:n))
                      ) + (1-w)*xold(i);
        end
        xnew
    end
    x = xnew;
end
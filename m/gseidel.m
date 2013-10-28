% ------------------------------------------------------------------------------
% FUNCTION:
%       gseidel
%
% PARAMS:
%       A - <nxn> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el sistema Ax=b mediante el metodo de Gauss-Seidel. No siempre 
%       hay convergencia.
% ------------------------------------------------------------------------------

function x = gseidel(A,b)
    [m n] = size(A);
    tol  = 0.0001;
    xold = ones(n,1)*10000;
    xnew = zeros(n,1);
    while ( norm(xnew-xold,"inf")/norm(xnew,"inf") > tol )
        xold = xnew;
        for i = 1:n
            xnew(i) = (1/A(i,i))*(
                        b(i)
                      - sum(A(i,1:i-1)*xnew(1:i-1))
                      - sum(A(i,i+1:n)*xnew(i+1:n))
                      );
        end
        xnew
    end
    x = xnew;
end
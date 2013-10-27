% ------------------------------------------------------------------------------
% FUNCTION:
%       qr_gschmidt
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       E - <mxn> numeric
%       U - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion QR de la matriz 'A'. Donde 'E = Q' y 'U = R'.
%       Se utiliza el metodo de Gram-Schmidt.
% ------------------------------------------------------------------------------

function [E U] = qr_gschmidt(A)
    [m n] = size(A);
    
    E = eye(m,n);
    U = zeros(n,n);
    for j = 1:n
        sum = A(:,j);
        for i = 1:j-1
            U(i,j) = dot(E(:,i),A(:,j));
            sum    = sum - U(i,j)*E(:,i);
        end
        U(j,j) = norm(sum);
        E(:,j) = sum/U(j,j);
    end
end
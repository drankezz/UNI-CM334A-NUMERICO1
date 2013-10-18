% ------------------------------------------------------------------------------
% FUNCTION:
%       doolitle
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       L - <nxn> numeric
%       U - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion LU de una matriz mediante el metodo de
%       Doolitle. Es posible que la matriz 'A' no tenga factorizacion LU por lo
%       que antes de iniciar el proceso de factorizacion se verifica utilizando
%       el lema 1.3 del libre 'Tecnicas de Calculo para Sistemas de Ecuaciones'
%       de Jose Luis de la Fuente O'connor.
% ------------------------------------------------------------------------------

function [L U] = doolittle(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' tenga factorización LU
    if !issquare(A)
        error("La matriz debe ser cuadrada.");
        return;
    end

    % Verifica que A admite factorizacion LU
    for i = 1:n
        if det(A(1:i,1:i)) == 0
            error("La matriz no tiene fact. LU. intenta pivotar la matriz\n");
            return;
        end
    end

    L = eye(n,n);
    U = eye(n,n);

    for k = 1:n
        for i = 1:k
            U(i,k) = A(i,k) - L(i,1:i-1)*U(1:i-1,k);

        end
        for i = k+1:n
            L(i,k) = (A(i,k) - L(i,1:k-1)*U(1:k-1,k))/U(k,k);
        end

    end
end

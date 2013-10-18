% ------------------------------------------------------------------------------
% FUNCTION:
%       croutl1u
%
% PARAMS:
%       A - <nxm> numeric
%
% RETURN:
%       L - <nxn> numeric
%       U - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion LU de Crout version L1U lo que quiere decir que
%       la diagonal de L son todos unos (1's). La matriz 'A' puede no tener
%       factorizacion LU por lo que se verifica utilizando el Lema 1.3 del libro
%       'Tecnicas de Calculo para Sistemas de Ecuaciones' de Jose Luis de la
%       fuente O'conor.
% ------------------------------------------------------------------------------

function [L U] = croutl1u(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' tenga factorización LU
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
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

    % Inicia la factorizacion LU
    for k = 1:n
        U(k,k:n)   = A(k,k:n) - L(k,1:k-1)*U(1:k-1,k:n);
        L(k+1:n,k) = (A(k+1:n,k) - L(k+1:n,1:k-1)*U(1:k-1,k))/U(k,k);
    end
end

% ------------------------------------------------------------------------------
% FUNCTION:
%       ldlt
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       L - <nxn> numeric
%       D - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion LDLT de una matriz simetrica 'A'. Donde 'L' es
%       una matriz triangular inferior unitaria (diagonal 1's) y 'D' es una
%       diagonal.
% ------------------------------------------------------------------------------

function [L D] = ldlt(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' tenga factorización LU
    if !issquare(A)
        error("La matriz debe ser cuadrada.");
        return;
    end

    % Comprueba que 'A' sea simetrica
    if !issymmetric(A)
        error("La matriz debe ser simetrica.");
        return;
    end

    L = eye(n);
    d = zeros(1,n);

    % Inicia la factorizacion LDLT
    for k = 1:n
        d(k) = A(k,k) - (L(k,1:k-1).^2)*d(1:k-1)';
        if d(k) == 0
            error('La matriz no es singular y no tiene factorizacion LDLT');
            return;
        end
        for i = k+1:n
            L(i,k) = (A(i,k) - sum(L(i,1:k-1).*L(k,1:k-1).*d(1:k-1)))/d(k);
        end
    end

    D = diag(d);
end

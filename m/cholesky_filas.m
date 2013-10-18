% ------------------------------------------------------------------------------
% FUNCTION:
%       cholesky_filas
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       G - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion de Cholesky G'G por filas de la matriz
%       simetrica 'A'. La matriz 'A' debe ser definida positiva por lo que se
%       comprobara antes de efectuar la factorizacion.
% ------------------------------------------------------------------------------

function G = cholesky_filas(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' sea cuadrada
    if !issquare(A)
        error("La matriz debe ser cuadrada.");
        return;
    end

    % Comprueba que 'A' sea simetrica
    if !issymmetric(A)
        error("La matriz debe ser simetrica.");
        return;
    end

    % Comprueba que 'A' sea definida positiva
    if isdefinite(A) <= 0
        error("La matriz debe ser definida positiva.");
        return;
    end

    G = eye(n);

    % Inicia el proceso de factorizacion Cholesky
    for i = 1:n
        G(i,i)     = sqrt(A(i,i) - sum(G(1:i-1,i).^2));
        for j = i+1:n   % De aqui el nombre de la variante 'por filas'
            G(i,j) = (A(i,j) - sum(G(1:i-1,i).*G(1:i-1,j)))/G(i,i);
        end
    end
end

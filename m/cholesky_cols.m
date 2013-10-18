% ------------------------------------------------------------------------------
% FUNCTION:
%       cholesky_cols
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       G - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion de Cholesky G'G por columnas de la matriz
%       simetrica 'A'. La matriz 'A' debe ser definida positiva por lo que se
%       comprobara antes de efectuar la factorizacion.
% ------------------------------------------------------------------------------

function G = cholesky_cols(A)
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
    for j = 1:n
        for i = 1:j-1   % De aqui el nombre de la variante 'por columnas'
            G(i,j) = (A(i,j) - sum(G(1:i-1,i).*G(1:i-1,j)))/G(i,i);
        end
        G(j,j) = sqrt(A(j,j) - sum(G(1:j-1,j).^2));
    end
end

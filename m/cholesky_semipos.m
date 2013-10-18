% ------------------------------------------------------------------------------
% FUNCTION:
%       cholesky_semipos
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       G - <nxn> numeric
%
% DESCRIPTION:
%       Factoriza la matriz 'A' por el metodo de Cholesky. La condicion es que
%       'A' puede ser semidefinida positiva. No se considera pivotacion.
% ------------------------------------------------------------------------------

function G = cholesky_semipos(A)
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

    for i = 1:n
        if A(i,i) > 0
            G(i,i) = sqrt(A(i,i) - sum(G(1:i-1,i).^2));
        end
        for j = i+1:n
            G(i,j) = (A(i,j) - sum(G(1:i-1,i).*G(1:i-1,j)))/G(i,i);
        end
    end
end

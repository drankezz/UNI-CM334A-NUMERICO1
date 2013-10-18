% ------------------------------------------------------------------------------
% FUNCTION:
%       sol_croutlu1_pparcial
%
% PARAMS:
%       A - <mxn> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%       F - <nxn> numerica
%
% DESCRIPTION:
%       Resuelve el sistema Ax=b mediante la factorizacion Crout LU1 de la
%       matriz 'A'. Se utiliza la pivotacion parcial para factorizar 'A' por lo
%       tambien se devuelve la matriz de pivotacion 'F'.
% ------------------------------------------------------------------------------

function [x F] = sol_croutlu1_pparcial(A, b)
    [m n] = size(A);

    % Comprueba que 'A' sea cuadrada
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

    % Comprueba que 'b' tiene el mismo numero de filas que columnas de 'A'
    if !iscolumn(b) || length(b) != n
        error("Las dimensiones son inconsistentes en el segundo argumento.");
    end

    [L U F] = croutlu1_pparcial(A);
    b = F*b;
    d = res_triang_inf(L,b);
    x = res_triang_sup(U,d);
end

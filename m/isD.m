% ------------------------------------------------------------------------------
% FUNCTION:
%       isD
%
% PARAMS:
%       A - <nxm> numeric
%
% RETURN:
%       r - numeric
%
% DESCRIPTION:
%       Verifica que 'A' sea una matriz diagonal en cuyo caso la funcion
%       devuelve 'r' = 1 si no devuelve 'r' = 0.
% ------------------------------------------------------------------------------

function r = isD(A)
    % Obtiene las dimensiones de 'A'
    [m n] = size(A);

    % Verifica que 'A' sea una matriz cuadrada
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
        return;
    end

    % Recorre todos los elementos de la matriz
    for f = 1:n
        for c = 1:n
            if f!=c && A(f,c) != 0
                r = 0;  % 'A' no es triangular inferior
                return;
            end
        end
    end
    r = 1;  % 'A' es una matriz triangular inferior
end

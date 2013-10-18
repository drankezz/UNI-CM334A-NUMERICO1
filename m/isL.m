% ------------------------------------------------------------------------------
% FUNCTION:
%       isL
%
% PARAMS:
%       A - <nxm> numeric
%
% RETURN:
%       r - numeric
%
% DESCRIPTION:
%       Verifica que 'A' sea una matriz cuadrada triangular Inferior en cuyo
%       caso la funcion devuelve 'r' = 1 si no devuelve 'r' = 0.
% ------------------------------------------------------------------------------

function r = isL(A)
    % Obtiene las dimensiones de 'A'
    [m n] = size(A);

    % Verifica que 'A' sea una matriz cuadrada
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
        return;
    end

    % Recorre los elementos superiores de la matriz 'A'
    for f = 1:n-1
        for c = f+1:n
            if A(f,c)
                r = 0;  % 'A' no es triangular inferior
                return;
            end
        end
    end
    r = 1;  % 'A' es una matriz triangular inferior
end

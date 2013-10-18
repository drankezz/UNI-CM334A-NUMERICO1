% ------------------------------------------------------------------------------
% FUNCTION:
%       isU
%
% PARAMS:
%       A - <nxm> numeric
%
% RETURN:
%       r - numeric
%
% DESCRIPTION:
%       Verifica que 'A' sea una matriz cuadrada triangular superior en cuyo
%       caso la funcion devuelve 'r' = 1 si no devuelve 'r' = 0.
% ------------------------------------------------------------------------------

function r = isU(A)
    % Obtiene las dimensiones de 'A'
    [m n] = size(A);

    % Verifica que 'A' sea una matriz cuadrada
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
        return;
    end

    % Recorre los elementos inferiores de la matriz 'A'
    for f = 2:n
        for c = 1:f-1
            if abs(A(f,c)) != 0 % Es posible que aparezca un -0
                r = 0;  % 'A' no es triangular superior
                return;
            end
        end
    end
    r = 1;  % 'A' es una matriz triangular superior
end

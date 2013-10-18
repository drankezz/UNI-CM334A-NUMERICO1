% ------------------------------------------------------------------------------
% FUNCTION:
%       mat2dmax
%
% PARAMS:
%       A - <nxm> numeric
%
% RETURN:
%       x - <1x3> numeric
%
% DESCRIPTION:
%       Encuentra el maximo valor dentro de la matriz 'A' (de dimension dos como
%       mucho) y lo devuelve, ademas del indice de su posicion dentro de la
%       matriz.
% ------------------------------------------------------------------------------

function [v f c] = mat2dmax(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    v = A(1,1);
    f = 1;
    c = 1;

    for i = 1:m
        for j = 1:n
            if(A(i,j) > v)
                v = A(i,j);
                f = i;
                c = j;
            end
        end
    end
end

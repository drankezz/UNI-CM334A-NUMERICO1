% ------------------------------------------------------------------------------
% FUNCTION:
%       res_triang_inf
%
% PARAMS:
%       A - <nxm> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el sistema Ax = b donde 'A' es una matriz cuadrada triangular
%       inferior.
% ------------------------------------------------------------------------------

function x = res_triang_inf(A, b)
    % Obtiene informacion de 'A'
    [m n] = size(A);

    % Comprueba que A sea cuadrada y triangular inferior
    if !issquare(A) || !isL(A)
        error("El argumento #1 no es una matriz cuadrada triangular inferior.");
        return;
    end

    % Resuelve el sistema de arriba hacia abajo
    x = zeros(n,1);
    for f = 1:n
        x(f) = (b(f) - A(f,1:f-1)*x(1:f-1))/A(f,f);
    end
end

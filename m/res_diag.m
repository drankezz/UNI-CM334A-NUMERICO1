% ------------------------------------------------------------------------------
% FUNCTION:
%       res_diag
%
% PARAMS:
%       A - <nxm> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el sistema Ax = b donde 'A' es una matriz diagonal.
% ------------------------------------------------------------------------------

function x = res_diag(A, b)
    % Obtiene informacion de 'A'
    [m n] = size(A);

    % Comprueba que A sea cuadrada y triangular inferior
    if !issquare(A) || !isD(A)
        error("El argumento #1 no es una matriz cuadrada diagonal.");
        return;
    end

    % Resuelve el sistema de arriba hacia abajo
    x = zeros(n,1);
    for f = 1:n
        x(f) = b(f)/A(f,f);
    end
end

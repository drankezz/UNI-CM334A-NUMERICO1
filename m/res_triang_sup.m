% ------------------------------------------------------------------------------
% FUNCTION:
%       res_triang_sup
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
%       superior.
% ------------------------------------------------------------------------------

function x = res_triang_sup(A, b)

    % Obtiene informacion de 'A'
    [m n] = size(A);

    % Comprueba que A sea cuadrada y triangular superior
    if !issquare(A) || !isU(A)
        error("El argumento #1 no es una matriz cuadrada triangular superior.");
        return;
    end

    % Resulve el sistema de abajo hacia arriba
    x = zeros(n,1);
    for f = n:-1:1
        x(f) = (b(f) - A(f,f+1:n)*x(f+1:n))/A(f,f);
    end
end

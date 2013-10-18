% ------------------------------------------------------------------------------
% FUNCTION:
%       gjordan_pparcial.m
%
% PARAMS:
%       A - <nxm> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el sistema Ax = b mediante la tecnica de Gauss-Jordan. El
%       algoritmo de eliminación emplea la tecnica de pivoteo parcial.
% ------------------------------------------------------------------------------

function x = gjordan_pparcial(A, b)
    % Obtiene las dimensiones de 'A'
    [m n] = size(A);

    % Comprueba que 'A' sea cuadrada
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
        return;
    end

    % Comprueba que 'b' tiene el mismo numero de filas que columnas de 'A'
    if !iscolumn(b) || length(b) != n
        error("Las dimensiones son inconsistentes en el segundo argumento.");
    end

    temp = eye(n);  % Matriz inicial para la diagonalizacion Gauss-Jordan

    % Crea la matrix extendida E
    E = [A b];

    % Inicia la eliminación por filas
    for i = 1:n

        % Inicio del proceso de pivotacion parcial
        [maxv maxi] = max(E(i:n,i));
        maxi        = i + maxi - 1;
        ftemp       = E(maxi,:);
        E(maxi,:)   = E(i,:);
        E(i,:)      = ftemp;
        % Fin del proceso de pivotacion parcial

        % f_j <- f_j - r*f_i
        for j = 1:n
            r = E(j,i)/E(i,i);
            if j == i
                r = 0;
            end
            E(j,:) = E(j,:) - r*E(i,:);
        end
        temp(i,i) = E(i,i);
    end
    % Termina la eliminacion por filas (Gauss)

    A = temp;   % Ahora 'A' es una matriz diagonal
    b = E(:,n+1);

    x = res_diag(A,b);
end

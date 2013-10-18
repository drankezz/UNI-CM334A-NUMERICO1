% ------------------------------------------------------------------------------
% FUNCTION:
%       gauss_pparcial
%
% PARAMS:
%       A - <nxm> numeric
%       b - <nx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el sistema Ax = b mediante la tecnica de eliminación de Gauss.
%       El algoritmo de eliminación de Gauss utilizado en este programa utiliza
%       la técnica de pivoteo parcial.
% ------------------------------------------------------------------------------

function x = gauss_pparcial(A, b)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' sea cuadrada
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
        return;
    end

    % Comprueba que 'b' tiene el mismo numero de filas que columnas de 'A'
    if !iscolumn(b) || length(b) != n
        error("Las dimensiones son inconsistentes en el segundo argumento.");
    end

    E = [A b];          % Matriz extendida

    % Inicia la eliminación por filas (Gauss) con pivotacion parcial
    for i = 1:n

        % Inicio del proceso de pivotacion parcial
        [maxv maxi] = max(E(i:n,i));
        maxi        = i + maxi - 1;
        temp        = E(maxi,:);
        E(maxi,:)   = E(i,:);
        E(i,:)      = temp;
        % Fin del proceso de pivotacion parcial

        % f_j <- f_j - r*f_j
        for j = i+1:n
            r = E(j,i)/E(i,i);
            E(j,:) = E(j,:) - r*E(i,:);
            E(j,i) = 0; % Para asegurar que no haya error de precision
        end

    end
    % Termina la eliminacion por filas (Gauss)

    A = E(:,1:n);   % Ahora 'A' es una matriz triangular superior
    b = E(:,n+1);
    x = res_triang_sup(A, b);
end

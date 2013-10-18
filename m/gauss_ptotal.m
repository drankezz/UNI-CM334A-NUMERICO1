% ------------------------------------------------------------------------------
% FUNCTION:
%       gauss_ptotal
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
%       la técnica de pivoteo total. Dado que cada columna contiene los
%       coeficientes de una variable, se hara seguimiento de las columnas para
%       mostrar la solucion en el orden correcto.
% ------------------------------------------------------------------------------

function x = gauss_ptotal(A, b)
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

    E  = [A b];         % Matriz extendida
    cc = 1:n;           % Guarda los intercambios entre columnas

    % Inicia la eliminación por filas (Gauss) con pivotacion parcial
    for i = 1:n
        % Inicio del proceso de pivotacion total
        [maxv maxf maxc] = mat2dmax(E(i:n,i:n));
        maxf = maxf + i - 1;
        maxc = maxc + i - 1;

        % f_i <-> f_maxf (Intercambio de filas)
        tempf     = E(maxf,:);
        E(maxf,:) = E(i,:);
        E(i,:)    = tempf;

        % c_i <-> c_maxc (Intercambio de columnas)
        tempc     = E(:,maxc);
        E(:,maxc) = E(:,i);
        E(:,i)    = tempc;

        % Seguimiento de las columas originales
        cc([i maxc]) = cc([maxc i]);
        % Fin del proceso de pivotacion total

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

    % Reordenando las filas de la solucion en 'temp'
    j = 1;
    temp = zeros(n,1);
    for i = cc
        temp(i) = x(j);
        j = j+1;
    end

    x = temp;
end

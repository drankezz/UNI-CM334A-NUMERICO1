% ------------------------------------------------------------------------------
% FUNCTION:
%       gjordan_ptotal.m
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
%       algoritmo de eliminación emplea la tecnica de pivoteo total. Dado que
%       cada columna contiene los coeficientes de una variable, se hara
%       seguimiento de las columnas para mostrar la solucion en el orden
%       correcto.
% ------------------------------------------------------------------------------

function x = gjordan_ptotal(A, b)
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
    cc = 1:n;           % Guarda los intercambios entre columnas

    % Inicia la eliminación por filas
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
    x = res_diag(A, b);

    % Reordenando las filas de la solucion en 'temp'
    j = 1;
    temp = zeros(n,1);
    for i = cc
        temp(i) = x(j);
        j = j+1;
    end

    x = temp;
end

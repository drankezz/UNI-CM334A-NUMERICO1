% ------------------------------------------------------------------------------
% FUNCTION:
%       croutlu1_pparcial
%
% PARAMS:
%       A - <nxm> numeric
%
% RETURN:
%       L - <nxn> numeric
%       U - <nxn> numeric
%       F - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion LU de Crout version LU1 lo que quiere decir que
%       la diagonal de U son todos unos (1's). Se utiliza pivotacion parcial con
%       lo cual se devuelve la  matriz de Permutacion de Filas 'F'. La matriz
%       'A' puede no tener factorizacion LU por lo que se verifica utilizando el
%       Lema 1.3 del libro 'Tecnicas de Calculo para Sistemas de Ecuaciones' de
%       Jose Luis de la fuente O'conor.
% ------------------------------------------------------------------------------

function [L U F] = croutlu1_pparcial(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' sea cuadrada
    if !issquare(A)
        error("El primer argumento debe ser una matriz cuadrada.");
        return;
    end

    % Verifica que A admite factorizacion LU
    for i = 1:n
        if det(A(1:i,1:i)) == 0
            error("La matriz no tiene fact. LU. intenta pivotar la matriz\n");
            return;
        end
    end

    L = zeros(n,n);
    U = eye(n,n);
    F = eye(n,n);

    % Inicia la factorizacion LU
    for k = 1:n
        L(k:n,k)   = A(k:n,k) - L(k:n,1:k-1)*U(1:k-1,k);

        % Inicio del proceso de pivotacion parcial
        [maxv maxi] = max(L(k:n,k));
        maxi        = k + maxi - 1;
        ftemp       = A(maxi,:); Ftemp     = F(maxi,:); Ltemp     = L(maxi,:);
        A(maxi,:)   = A(k,:)   ; F(maxi,:) = F(k,:)   ; L(maxi,:) = L(k,:);
        A(k,:)      = ftemp    ; F(k,:)    = Ftemp    ; L(k,:)    = Ltemp;
        % Fin del proceso de pivotacion parcial

        U(k,k+1:n) = (A(k,k+1:n) - L(k,1:k-1)*U(1:k-1,k+1:n))/L(k,k);
    end
end

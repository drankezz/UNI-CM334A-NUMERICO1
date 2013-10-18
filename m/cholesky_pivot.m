% ------------------------------------------------------------------------------
% FUNCTION:
%       cholesky_pivot
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       G - <nxn> numeric
%       F - <nxn> numeric
%
% DESCRIPTION:
%       Realiza la factorizacion de Cholesky G'G por filas de la matriz
%       simetrica 'A'. La matriz 'A' debe ser definida positiva por lo que se
%       comprobara antes de efectuar la factorizacion. Se considera pivotaciones
%       Por lo que tambien se devuelve la matriz 'F' de permutaciones.
% ------------------------------------------------------------------------------

function [G F] = cholesky_pivot(A)
    [m n] = size(A);    % dimensiones de A (filas, columnas)

    % Comprueba que 'A' sea cuadrada
    if !issquare(A)
        error("La matriz debe ser cuadrada.");
        return;
    end

    % Comprueba que 'A' sea simetrica
    if !issymmetric(A)
        error("La matriz debe ser simetrica.");
        return;
    end

    % Comprueba que 'A' sea definida positiva
    if isdefinite(A) <= 0
        error("La matriz debe ser definida positiva.");
        return;
    end

    G = zeros(n);
    F = eye(n);

    % Inicia el proceso de factorizacion Cholesky
    for i = 1:n
        [maxv maxi] = max(diag(A(i:n,i:n)));
        maxi = maxi + i - 1;
        if A(maxi, maxi) > 0
            % Inicio de pivotacion
            % Intercambia las fila 'i' y 'maxi'
            ftemp     = A(maxi,:); Gtemp     = G(maxi,:); Ftemp     = F(maxi,:);
            A(maxi,:) = A(i,:)   ; G(maxi,:) = G(i,:)   ; F(maxi,:) = F(i,:)   ;
            A(i,:)    = ftemp    ; G(i,:)    = Gtemp    ; F(i,:)    = Ftemp    ;

            % Intercambia las columna 'i' y 'maxi'
            ctemp     = A(:,maxi); Gtemp     = G(:,maxi);
            A(:,maxi) = A(:,i)   ; G(:,maxi) = G(:,i)   ;
            A(:,i)    = ctemp    ; G(:,i)    = Gtemp    ;
            % Fin de pivotacion

            % Inicia la factorizacion de Cholesky
            G(i,i) = sqrt(A(i,i) - sum(G(1:i-1,i).^2));
            for j = i+1:n
                G(i,j) = (A(i,j) - sum(G(1:i-1,i).*G(1:i-1,j)))/G(i,i);
            end
        end
    end

    % Se puede comprobar mediante
    % F*G'*G*F
    % O sea, pivotar G'*G intercambiando filas y columnas
end

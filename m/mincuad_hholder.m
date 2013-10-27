% ------------------------------------------------------------------------------
% FUNCTION:
%       mincuad_hholder
%
% PARAMS:
%       A - <mxn> numeric
%       b - <mx1> numeric
%
% RETURN:
%       x - <nx1> numeric
%
% DESCRIPTION:
%       Resuelve el problema de minimos cuadrados Ax=b por el metodo de
%       Householder.
% ------------------------------------------------------------------------------

function x = mincuad_hholder(A,b)
    [m n] = size(A);
    [Q R] = qr_hholder(A);
    c = Q*b;
    c = c(1:n,1);   % Ahora c tiene dimensiones <nx1>
    
    %% Obteniendo la matriz R1
    R1 = zeros(n,n);
    for i = 1:n
        for j = 1:n
            if (j >= i)
                R1(i,j) = R(i,j);
            end
        end
    end
    
    x = res_triang_sup(R1,c);
end
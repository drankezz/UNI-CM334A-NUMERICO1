% ------------------------------------------------------------------------------
% FUNCTION:
%       qr_givens
%
% PARAMS:
%       A - <mxn> numeric
%       
% RETURN:
%       Q - <mxm> numeric
%       R - <nxn> numeric
%
% DESCRIPTION:
%       Realiza el proceso QR a la matriz 'A' por el metodo de Givens. 'Q' es 
%       una matriz ortogonal. 'R' es una matriz triangular superior. Se ha 
%       priorizado la claridad del codigo sobre la eficienci del mismo.
% ------------------------------------------------------------------------------

function [Q R] = qr_givens(A)
    [m n] = size(A);
    Q = eye(m,m);
    for j = 1:n
        for i = j+1:m
            %% Hacer nulo el elemento (i,j)
            G = eye(m,m); %% Matriz de rotacion de Givens para anular A(i,j)
            if (A(j,j) != 0)
                theta = atan(A(i,j)/A(j,j));
                c = cos(theta);
                s = sin(theta);
                G(i,i) =  c;
                G(j,j) =  c;
                G(i,j) = -s;
                G(j,i) =  s;
                A = G*A; %% Aplicando la transf. de Givens a 'A'
                Q = G*Q;
            end
        end
    end
    R = A;
end
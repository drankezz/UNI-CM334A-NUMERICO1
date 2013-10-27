% ------------------------------------------------------------------------------
% FUNCTION:
%       qr_hholder
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       Q - <mxm> numeric
%       R - <nxn> numeric
%
% DESCRIPTION:
%       Hace la triangularizacion ortogonal a la matriz 'A'. 'Q' es una matriz 
%       ortogonal y 'R' es una matriz triangular superior. Se utiliza el metodo 
%       de Householder. Para comprobar: Q*A = R
% ------------------------------------------------------------------------------

function [Q R] = qr_hholder(A)
    [m n] = size(A);
    
    Q = eye(m,m);
    
    w = zeros(m,n);
    for j = 1:n
    
        %% Proceso para hallar el vector w de householder de A(:,j)
        w(j:m,j) = A(j:m,j);
        s        = norm(w(:,j));
        r        = 1/sqrt(2*s*(s+abs(A(j,j))));
        w(j,j)   = A(j,j) + s*sign(A(j,j));
        w(:,j)   = r*w(:,j);
        
        %% Transformacion de householder para A(:,j)
        H = eye(m,m) - 2*w(:,j)*w(:,j)';
        
        %% Aplicacion de la transformacion de householder a A
        A = H*A;
        
        %% Actualizando Q
        Q = H*Q;
    end
    R = A;
end
% ------------------------------------------------------------------------------
% FUNCTION:
%       gschmidt
%
% PARAMS:
%       A - <mxn> numeric
%
% RETURN:
%       U - <mxn> numeric
%
% DESCRIPTION:
%       Realiza la ortonormalizacion de Gram-Schmidt en base a los vectores
%       columna de la matriz 'A'. Los vectores columna de 'e' son ortonormales 
%       entre si
% ------------------------------------------------------------------------------

function e = gschmidt(A)
    [m n] = size(A);
    
    e = eye(m,n);
    for j = 1:n
        sum = A(:,j);
        for i = 1:j-1
            sum = sum - dot(A(:,j),e(:,i))*e(:,i);
        end
        e(:,j) = sum/norm(sum);
    end
end
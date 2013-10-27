% ------------------------------------------------------------------------------
% FUNCTION:
%       mincuad_gschmidt
%
% PARAMS:
%       A - <mxn> numeric
%       b - <mx1> numeric
%
% RETURN:
%       x - <nx1>
%
% DESCRIPTION:
%       Resuelve el problema de minimos cuadrados mediante el metodo de 
%       Gram-Schmidt.
% ------------------------------------------------------------------------------

function x = mincuad_gschmidt(A,b)
    [E U] = qr_gschmidt(A);
    x = res_triang_sup(U, E'*b);
end
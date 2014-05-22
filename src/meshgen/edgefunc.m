function d = edgefunc( x, y )
% EDGEFUNC
%
% d = EDGEFUNC(x,y) computes the `signed distance' from each node located in
% (x,y) to a given curve.

    d = abs( sum( sign( y - line(x) ) ) );

end

function y = line( x )
%Y = LINE( X )
%
% Line describing the boundary.  (Usually y = 0).
%
    % y = 1.3*x;
    y = zeros( size(x) );

end

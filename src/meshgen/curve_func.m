function z = curve_func(xvec)
%CURVE_FUNC.  Zero set of intersection curve.
% z = CURVE_FUNC( xvec ) function defining the curve.  
% The set of all points satisfying curve_func(xpts) = 0 is the curve of
% interest.

    x = xvec(:,1);
    y = xvec(:,2);

    z = y - (-0.5 + x.^2);

end

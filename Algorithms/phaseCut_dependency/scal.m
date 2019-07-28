function r = scal (v1,v2)

% Scalar product of v1 and v2.

r = sum(conj(v1(:)).*v2(:)) ;
% -------------------------------------
function rect_vect = rectify(in_vect)
% function rect_vect = rectify(in_vect)
%
% Support function to rectify (eliminate negative values from) a vector
% ----------------------------------------
rect_vect = in_vect;
rect_vect(rect_vect < 0) = 0; 
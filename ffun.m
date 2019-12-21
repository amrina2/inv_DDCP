%%%%%%
%%%%% File name: ffun.m 
%%%%% Computing Artifact
%%%%% Author: Amrina Ferdous
%%%%% Purpose: This is the file where researchers need to formulate the
%%%%% 'f' function so that we could use it in the main.m file for
%%%%% calculating the frontier of the geophysics ill-posed nonlinear 
%%%%% inverse problem. Here 'f' represents the integrand in the forward
%%%%% operator.
%%%%%%%%%%%%%%%%%%%%%%%
%
function f = ffun(x,w,z);
% x= scalar
% w= vector with dimension m
% x are in n data point and it should be in a column
% w is discreatize with m point and should be in a row
%%%%%%%
 H = 10; % Deepth between the surface and the subsurface 
 m= length(w)
 n= length(x)
 x= repmat(x,1,m)
 w= repmat(w,n,1)
 z=repmat(z,n,1)
 f=log(((x-w).^2 + H.^2)./ ((x-w).^2 + (H-z).^2));
end
%

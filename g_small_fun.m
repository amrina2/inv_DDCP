function g = g_small_fun(x,w,z)
% x= scalar
% w= vector with dim m
% x are in n data point & it should be in a colums
% w is discreatize with m point & should be in a row
 H = 10;
 m= length(w);
 n= length(x);
 x1= repmat(x,1,m);
 w1= repmat(w,n,1);
 z1=repmat(z,n,1);
 f=log(((x1-w1).^2 + H.^2)./ ((x1-w1).^2 + (H-z1).^2));
 g=trapz(w,f,2);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


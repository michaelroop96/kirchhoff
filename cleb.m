function [a,b,c]=cleb(N)
a = rand(N,1);
c22 = rand(1);
c33 = rand(1);
c11 = a(2)*a(3)/(a(2)-a(3))*(c22/a(3)-c33/a(2)-(c22-c33)/a(1));
b11 = rand(1);
b22 = b11;
b33 = b11;
b = diag([b11,b22,b33]);
c = diag([c11,c22,c33]);
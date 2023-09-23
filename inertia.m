function [M1,M2]=inertia(W,theta,a,B,C)
m=[W(3,2);W(1,3);W(2,1)];
p=[theta(3,2);theta(1,3);theta(2,1)];
w=zeros(3,1);
Bp=B*p;
for i=1:3
    w(i)=a(i)*m(i)+Bp(i);
end
u=B*m+C*p;
M1 = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
M2 = [0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];



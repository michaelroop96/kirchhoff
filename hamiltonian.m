function H=hamiltonian(a,b,c,W,theta)
N = length(a);
H1 = 0;
H2 = 0;
H3 = 0;
m=[W(3,2);W(1,3);W(2,1)];
p=[theta(3,2);theta(1,3);theta(2,1)];
for i=1:N
    H1=H1+a(i)*m(i)^2;
end
for k=1:N
    for j=1:N
        H2=H2+b(k,j)*(p(k)*m(j)+m(k)*p(j));
        H3=H3+c(k,j)*p(k)*p(j);
    end
end
H=1/2*(H1+H2+H3);
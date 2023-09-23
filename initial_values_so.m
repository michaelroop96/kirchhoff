function X0 = initial_values_so(NMAT)

W0 = rand(NMAT,NMAT);
%W0 = rand(NMAT,NMAT);
W0 = (W0-W0')/2;
%W0 = W0-trace(W0)*eye(NMAT)/NMAT;
X0 = W0/norm(W0);

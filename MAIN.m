clear all
close all
clc
h = 1e-1;%Step of time discretization
NMAT = 3;%Dimension of matrices     
itermax = 10000; %Number of iterations
%[a,b,c]=generic(NMAT);
%[a,b,c]=kir(NMAT);
%[a,b,c]=cleb(NMAT);
[a,b,c]=LSK(NMAT);
B = (b+b')/2;
C = (c+c')/2;
W0 = initial_values_so(NMAT);
theta0 = initial_values_so(NMAT);
lambda0 = sort(imag(eig(theta0)));
k=1;
cas0 = trace(W0*theta0^k);
H0 = hamiltonian(a,b,c,W0,theta0);
dlambda = zeros(1,itermax);
cas = zeros(1,itermax);
VAR_HAM = zeros(1,itermax);
W = W0;
theta = theta0;
COMPW = zeros(itermax,NMAT,NMAT);
COMPTh = zeros(itermax,NMAT,NMAT);
for i=1:1:NMAT
        for j=1:1:NMAT
            COMPW(1,i,j) = W(i,j);
            COMPTh(1,i,j) = theta(i,j);
        end
end
for iter=1:itermax-1    
    VAR_HAM(iter)=real(hamiltonian(a,b,c,W,theta)-H0);
    dsp = sort(imag(eig(theta)))-lambda0;
   %dsp = sort(imag(eig(W)))-lambda0;
    dlambda(iter)=dsp(1);
    cas(iter) = real(trace(W*theta^k)-cas0);
    [tildeW,tildeTheta,nev1,nev2] = fixed_point_iteration(W,theta,NMAT,h,a,B,C);
%     tildeW = tildeW-trace(tildeW)*eye(NMAT)/NMAT;
%     tildeTheta = tildeTheta-trace(tildeTheta)*eye(NMAT)/NMAT;
    [M1,M2]=inertia(tildeW,tildeTheta,a,B,C);
    W = W+h*LieBracket(tildeW,M1)+h*LieBracket(tildeTheta,M2);
    theta = theta+h*LieBracket(tildeTheta,M1);
    for i=1:1:NMAT
        for j=1:1:NMAT
            COMPW(iter+1,i,j) = W(i,j);
            COMPTh(iter+1,i,j) = theta(i,j);
        end
    end
%     theta = (eye(NMAT)-h/2*M1)*tildeTheta/(eye(NMAT)-h/2*M1);
%     theta = theta-trace(theta)*eye(NMAT)/NMAT;

end

VAR_HAM(itermax)=real(hamiltonian(a,b,c,W,theta)-H0);
dsp = sort(imag(eig(theta)))-lambda0;
%dsp = sort(imag(eig(W)))-lambda0;
dlambda(itermax)=dsp(1);
cas(itermax) = real(trace(W*theta^k)-cas0);
plot_graphics(dlambda,cas,VAR_HAM,COMPW,COMPTh,itermax,h);
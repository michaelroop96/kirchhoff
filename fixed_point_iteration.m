function [tildeW,tildeTheta,nev1,nev2]=fixed_point_iteration(W,theta,NMAT,h,a,B,C)
%tol = 10^(-3);
W0 = zeros(NMAT,NMAT);
theta0 = zeros(NMAT,NMAT);
for i=1:30
    
    [M1,M2]=inertia(W0,theta0,a,B,C);
    W1 = W0-(-W+W0-h/2*LieBracket(W0,M1)-h/2*LieBracket(theta0,M2)-h^2/4*(M1*W0*M1+M2*theta0*M1+M1*theta0*M2));
    theta1 = theta0-(-theta+theta0-h/2*LieBracket(theta0,M1)-h^2/4*M1*theta0*M1);
%     W1 = W0-(-W+W0-h/2*LieBracket(W0,M1)-h/2*LieBracket(theta0,M2)-h^2/4*(M1*W0*M1+M1*M2*theta0+M1*theta0*M2));
%     theta1 = theta0-(-theta+(eye(NMAT)+h/2*M1)*theta0/(eye(NMAT)+h/2*M1));
    W0 = W1;
    theta0 = theta1;
   
    
end

tildeW = (W1-W1')/2;
tildeTheta = (theta1-theta1')/2;

[M1,M2]=inertia(tildeW,tildeTheta,a,B,C);
nev1 = norm(-W+tildeW-h/2*LieBracket(tildeW,M1)-h/2*LieBracket(tildeTheta,M2)-h^2/4*(M1*tildeW*M1+M2*tildeTheta*M1+M1*tildeTheta*M2));
nev2 = norm(-theta+tildeTheta-h/2*LieBracket(tildeTheta,M1)-h^2/4*M1*tildeTheta*M1);






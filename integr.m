function [nip2,xip2,w2,Nextr] = integr()
% Quadrature + extrapolation (nelvert=6 ONLY)

% 4-pt triangle rule (same as your code)
al=.6; be=.2;
nip2=4;
xip2=[al be be 1/3;
      be al be 1/3];
w2  =[25,25,25,-27]/48;

mxi=[1  0 1;
     0  1 1;
     0  0 1;
     .5 .5 1;
     0 .5 1;
     .5  0 1];

% build N and extrap operator
N = zeros(nip2,6);
for i2=1:nip2
    xi1=[xip2(:,i2); 1-sum(xip2(:,i2))];
    xi11=xi1(1); xi12=xi1(1)-.5;
    xi21=xi1(2); xi22=xi1(2)-.5;
    xi31=xi1(3); xi32=xi1(3)-.5;
    N(i2,:)=2*[xi11*xi12; xi21*xi22; xi31*xi32; 2*xi11*xi21; 2*xi21*xi31; 2*xi11*xi31]';
end
Np    = pinv(N);
hxi   = [xip2' ones(nip2,1)];
hxip  = pinv(hxi);
Nextr = Np*(eye(nip2)-hxi*hxip)+mxi*hxip;
end

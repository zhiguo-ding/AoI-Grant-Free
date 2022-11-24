close all
K=2;
p=[0.001: 0.00001:1];
p1=1/K;
M=10;

y=(1-p+p*p1).^(M-1).*p*p1 + (1-p+(M-1)*p*p1).*(1-p).^(M-2).*p*p1;
aoi = (2-y)./y;
[a1,a2] = min(aoi);


pp=K/M;
yz=(1-pp+pp*p1).^(M-1).*pp*p1 + (1-pp+(M-1)*pp*p1).*(1-pp).^(M-2).*pp*p1;
plot(p,  (2-y)./y,p,(2-yz)/yz*ones(1,length(p)))
ylim([0 5000])
%plot(p,y)


z=p(a2)/2
%z=2/M/2
(1-M*z)*(1-z)^(M-2)+(1-4*z-2*M*(M-3)*z^2)*(1-2*z)^(M-3)
[(sqrt(16+8*M*(M-3))-4)/4/M/(M-3) sqrt(2)/M/2 p(a2)/2 2/M/2]

exp(1)/(2+exp(1))
sqrt(2)*exp(sqrt(2)-1)/(exp(sqrt(2)/2)+1+sqrt(2)/2)

M=[1: 1: 100];
plot(M, (sqrt(4+2*M.*(M-1))-4)./M - 2./M)
% j=0;
% for z=0.1 :0.01:1
%     j= j+1;
%     ttt1=(1-M*z)*(1-z)^(M-2)+(1-4*z-2*M*(M-3)*z^2)*(1-2*z)^(M-3);
%     ttt=(1-z)^(M-1)-(M-1)*(1-z)^(M-2)*z+...
%         (1+2*(M-3)*z)*(1-2*z)^(M-2)-2*(M-2)*(z+(M-3)*z^2)*(1-2*z)^(M-3);
%     sx(j) = ttt1-ttt;
% end
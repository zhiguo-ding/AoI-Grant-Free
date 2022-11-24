clear all 


%M = 32;%number of users
%N = 8;%number of slot within one frame
K = 2; %number of power levels
T = 6;%duration of a slot
Mvec = [1 5 10 15 20  25 30];
Nvec = [1 5 10 15 20 25 30  ];
pjvec = [0.01:0.001:1];
P = 100;
R=0.5; eps = 2^R-1;
for iave = 1: length(pjvec) 
%M = Mvec(iave);    
M = 8;
%N = Nvec(iave);   
N = 1;%number of slot within one frame

J = 0;%number of frames 
num_frame = 50; %num_frame frames are used
h = complex(sqrt(0.5)*randn(M,N*num_frame),sqrt(0.5)*randn(M,N*num_frame));
h = abs(h).^2;

%the optimal value      
if M==1
    t_opt = fzero(@(t) (1-t/2)*exp(-t/2)+(1-t^2/2)*exp(-t), M*min(K,M)/M);
else
    t_opt = fzero(@(t) (1-t/2)*exp(-t/2)+(1-t^2/2)*exp(-t),[M*(sqrt(4+2*M*(M-1))-4)/M/(M-3), M*min(K,M)/M]);
end
if K==2
    p_tx = min(1,t_opt/M);
else
    p_tx = min(K/M,1);
end
p_tx = pjvec(iave);%fixed choice

S_j = [];
Y_j = [];
X_j = [];
t_succ = [];
t_last = 0;%the time for the last successful update
J=0; % record the number of successful updates
x_j = 0; 

%%%%%%%%%%%% decide those NOMA power levels
Pnoma = zeros(K,1);
Pnoma(K) = 2^R-1;

for k = K-1 :-1: 1
    Pnoma(k) = (2^R-1)*(1+sum(Pnoma(k+1:end)));
end
%%%%%%%%%%
zzz=0;
sum1=0; sum2=0;sum3=0; sum4=0; 

%analysis%%%%%%%%%%%%%%%%%%%%%%%%%
%find the P matrix
Pmatrix = zeros(M,M);
p1=1/K; 
%pjj
for j = 0: M-1
     p1=1/K;
     pjj=0;
     for m = 1: M-j
         cmm2 = factorial(M-j)/factorial(m)/factorial(M-j-m);
         sumtemp1 = 0;
         for k = 1 : K
             sumtemp1 = sumtemp1 + m*p1*(1-k*p1)^(m-1);
         end
         pjj = pjj + cmm2*p_tx^m*(1-p_tx)^(M-j-m)*sumtemp1;
     end
     Pmatrix(j+1,j+1) = 1-pjj;
end

%p_jj+1
for j = 0 : M-2
     pjj1 =0;
     for m = 2: M-j
         cmm2 = factorial(M-j)/factorial(m)/factorial(M-j-m);
         sumtemp1 = 0;
         for k = 1 : K
             sumtemp2 = 0;
             for ka1=k+1 : K
                 sumtemp2 = sumtemp2 + (m-1)*p1*(1-ka1*p1)^(m-2);
                 if ka1==K & (m-1)*p1*(1-ka1*p1)^(m-2)~=0
                     dfd=1;
                 end
             end

             sumtemp1 = sumtemp1 + (m/(M-j)*(m-1)+(M-j-m)/(M-j)*m)*p1...
                 *((1-k*p1)^(m-1) - sumtemp2);
         end
         pjj1 = pjj1 + cmm2*p_tx^m*(1-p_tx)^(M-j-m)*sumtemp1;
     end
     Pmatrix(j+1,j+2) = pjj1 + (M-j)*p_tx*(1-p_tx)^(M-j-1)*(M-j-1)/(M-j);
end

%for the rest, we can use the generic formula
for i = 2 : K %p_j j+i
     for j = 0 : M-i-1
         pjj1 =0;
         for m = i+1: M-j
             cmm2 = factorial(M-j)/factorial(m)/factorial(M-j-m);
             sumtemp1 = 0;
             for k1 = 1 : K
                 for k2 = k1+i-1 : K
                     sumtemp2 = 0;
                     for ka1=k2+1 : K
                         sumtemp2 = sumtemp2 + (m-i)*p1*(1-ka1*p1)^(m-i-1);
                     end
                     mpp = [m:-1:m-i+1];
                     cmk = factorial(k2-k1-1)/factorial(i-2)/factorial(k2-k1-i+1);
                     sumtemp1 = sumtemp1 + (M-j-i)/(M-j)*prod(mpp)*p1^i*cmk....
                         *((1-k2*p1)^(m-i) - sumtemp2);
                 end
             end
             pjj1 = pjj1 + cmm2*p_tx^m*(1-p_tx)^(M-j-m)*sumtemp1;
         end
         %a special case m=i
         m=i; 
             cmm2 = factorial(M-j)/factorial(m)/factorial(M-j-m);
             sumtemp1 = 0;
             for k1 = 1 : K
                 for k2 = k1+i-1 : K
                     sumtemp2 = 0;
                     mpp = [m:-1:m-i+1];
                     cmk = factorial(k2-k1-1)/factorial(i-2)/factorial(k2-k1-i+1);
                     sumtemp1 = sumtemp1 + (M-j-i)/(M-j)*prod(mpp)*p1^i*cmk....
                         *((1-k2*p1)^(m-i) - sumtemp2);
                 end
             end
             Pmatrix(j+1,j+i+1) = pjj1 + cmm2*p_tx^m*(1-p_tx)^(M-j-m)*sumtemp1;
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the p vector 
p_vector = zeros(M,1);
for j = 1 : M
    p_vector(j,1) = 1- sum(Pmatrix(j,:)) ;
end
    
s0 = zeros(M,1); s0(1)=1;
pfail = s0'*Pmatrix^N*ones(M,1);

EYj = T*N/(1-pfail);
EXj = 1/(1-pfail);
EXj2 = (1+pfail)/(1-pfail)^2;
ESj = 0; ESj2 = 0;
for n = 1 : N
    ESj = ESj + T*n*s0'*Pmatrix^(n-1)*p_vector/(1-pfail);
    ESj2 = ESj2 +T^2* n^2*s0'*Pmatrix^(n-1)*p_vector/(1-pfail);
end

EYj2 = T^2*N^2*EXj2 + 2*ESj2 - 2*(ESj)^2;
ESj1Y = ESj*EYj - ESj2 + (ESj)^2;

Aoi_ana(iave) = ESj1Y/EYj + EYj2/2/EYj;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% the benchmarking scheme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_j = [];
Y_j = [];
X_j = [];
t_succ = [];
t_last = 0;%the time for the last successful update
J=0; % record the number of successful updates
x_j = 0; 
p_tx =  0.05* ones(1,M);% 1./(M-[0:M-1]);%0.1* ones(1,M);%
p_tx =   1./(M-[0:M-1]);%
p_tx = pjvec(iave)* ones(1,M);
 
                     
         
%analysis
P_matrix = zeros(M,M);
for j = 0 : M-1
    for i = 0 : M-1
        if i==j
            P_matrix(j+1,i+1) = 1 - (M-j)*p_tx(j+1)*...
                exp(-eps/P)*(1-p_tx(j+1))^(M-j-1);
        elseif i == j+1
            P_matrix(j+1,i+1) = (M-j-1)*p_tx(j+1)*...
                exp(-eps/P)*(1-p_tx(j+1))^(M-j-1);
        end
    end
end

p_vector = zeros(M,1);
for j = 0 : M-1
    p_vector(j+1,1) = p_tx(j+1)*...
                exp(-eps/P)*(1-p_tx(j+1))^(M-j-1);
end
    
s0 = zeros(M,1); s0(1)=1;
pfail = s0'*P_matrix^N*ones(M,1);

EYj = N/(1-pfail);
EXj = 1/(1-pfail);
EXj2 = (1+pfail)/(1-pfail)^2;
ESj = 0; ESj2 = 0;
for n = 1 : N
    ESj = ESj + n*s0'*P_matrix^(n-1)*p_vector/(1-pfail);
    ESj2 = ESj2 + n^2*s0'*P_matrix^(n-1)*p_vector/(1-pfail);
end

EYj2 = N^2*EXj2 + 2*ESj2 - 2*(ESj)^2;
ESj1Y = ESj*EYj - ESj2 + (ESj)^2;

Aoi_ana_oma(iave) = T*ESj1Y/EYj + T*EYj2/2/EYj;
end

plot( pjvec,Aoi_ana_omaf, pjvec, Aoi_anaf) 




% ax2 = axes('Position',[.7 .7 .2 .2])
% box on;
%plot(Mvec,aoi_oma, Mvec,Aoi_ana_oma, '-s',Mvec,aoi, Mvec, Aoi_ana,'-o') 
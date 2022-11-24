clear all 


%M = 32;%number of users
%N = 8;%number of slot within one frame
K = 3; %number of power levels
T = 6;%duration of a slot
Mvec = [1 5 10  15 20 25   30 ];
Nvec = [1 5 10 15 20 25 30  ];
P = 100;
R=0.5; eps = 2^R-1;
for iave = 1: length(Mvec) 
M = Mvec(iave);    
%M = 8;
%N = Nvec(iave);   
N = 8;%number of slot within one frame

J = 0;%number of frames 
num_frame = 500000; %num_frame frames are used
h = complex(sqrt(0.5)*randn(M,N*num_frame),sqrt(0.5)*randn(M,N*num_frame));
h = abs(h).^2;

%the optimal value  
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
%p_tx = 0.05;%fixed choice

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
for i = 1 : num_frame %will use user 1 as the targged user
    succ = 0; 
    n = 1; %index of time slow within a frame
    m = 0; %number of successful devices
    while (succ ==0)      &   n<=N          &  m<M
        if n == 1%a new frame starts
           x_j = x_j + 1; 
           tempjj=0; 
        end
        %tagged user fails  within the frame .  not all user scceed
        h_current = h(1:M-m,(i-1)*N+n); %get the channels for M-m active users
        
        if M==1%for the special case M=1
            tx_level=K;
        else            
            tx_level = randsrc(M-m,1,[1: K]); %which level they want to use
        end
        tx_prob = binornd(1,p_tx(1)*ones(M-m,1)); %1 means tx: willness to send
        
        tx_activex = tx_level.*tx_prob; % activie users' level
        
        for k =1 : K %K power levels
            if M==1
                k=K;
            end
            tx_outage = (Pnoma(k)./h_current(:,1)<=P);
            tx_active = tx_activex.*tx_outage;
            setk = find(~(tx_active-k)); %the set contains users choose power level k
            if length(setk)==0 %no one chooses this level
                %nothing to do, go to the next level
            elseif length(setk)>1 %collision leads to SIC failure
                break;
            else %only one user chooses the level             

                        user_index = setk;
                        if user_index == 1 %the tagged user succeeds
                            succ = 1;
                            J = J +1; %one more successful updates
                            S_j = [S_j n];%record the access delay
                            Y_j = [Y_j (i-1)*N+n-t_last]; %update Y_j
                            t_last = (i-1)*N+n;
                            t_succ = [t_succ (i-1)*N+n];
                            X_j = [X_j x_j];
                            x_j = 0;%for each update, x_j is back to 0   
                            
                            break; % exit the loop with K and go to the next frame
                        else %a user other than the tagged user succeeds
                            
                            m = m +1;% add one active user
                        end                       

                    
            end 
            if k==K %for the special case M=1
                break;
            end
        end
        n=n+1;
    end

end
 %sum3/num_frame -sum4/num_frame 
 
 
 %probabily for a single user (not 1) succeeds - probabily for two users (any) succeed
Q_j = [];
for i = 1 : J-1
    Q_j = [Q_j (S_j(i)*Y_j(i+1) + Y_j(i+1)^2/2)*T^2 ];
end
aoi(iave) = sum(Q_j)/sum(Y_j(2:end)*T);

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
         for k = 1 : K-1
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
             for k1 = 1 : K-i
                 for k2 = k1+i-1 : K-1
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
             for k1 = 1 : K-i+1
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

for i = 1 : num_frame %will use user 1 as the targged user
    if i == 36
        dfd=1;
    end
    succ = 0; 
    n = 1; %index of time slow within a frame
    m = 0; %number of successful devices
    while (succ ==0)      &   n<=N          &  m<M
       %tagged user fails  within the frame .  not all user scceed
       if n == 1%a new frame starts
           x_j = x_j + 1; 
           tempjj=0; 
       end
        h_current = h(1:M-m,(i-1)*N+n); %get the channels for M-m active users
        tx_status = binornd(1,p_tx(m+1)*ones(M-m,1)); %1 means tx
        if sum(tx_status) == 1 % a single user decides to transmit
            temp1 =  (h_current>=eps/P).*tx_status;
            if sum(temp1)==1 %that user succeeds 
                user_index = find(temp1);
                if user_index == 1 %the tagged user succeeds
                    succ = 1;
                    J = J +1; %one more successful updates
                    S_j = [S_j n];%record the access delay
                    Y_j = [Y_j (i-1)*N+n-t_last]; %update Y_j
                    t_last = (i-1)*N+n;
                    t_succ = [t_succ (i-1)*N+n];
                    X_j = [X_j x_j];
                    x_j = 0;%for each update, x_j is back to 0   
                    tempjj = 1; 
                else %a user other than the tagged user succeeds
                    m = m +1;% add one active user
                end
                              
            end
        end
        n = n + 1; % move to the next time slot   
    end    

end

Q_j = [];
for i = 1 : J-1
    Q_j = [Q_j (S_j(i)*Y_j(i+1) + Y_j(i+1)^2/2)*T^2 ];
end
aoi_oma(iave) = sum(Q_j)/sum(Y_j(2:end)*T);
                     
         
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

%plot(Nvec,aoi_oma, Nvec,Aoi_ana_oma, '-s',Nvec,aoi, Nvec, Aoi_ana,'-o') 
plot(Mvec,aoi_oma, Mvec,Aoi_ana_oma, '-s',Mvec,aoi, Mvec, Aoi_ana,'-o') 
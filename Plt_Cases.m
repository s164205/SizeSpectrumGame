LSD = 1/17; %log-sate-difference
M = exp((1:200)*LSD)/exp(LSD);

dt = .03;
T = 1:dt:50;

offset = round(log(400)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
omega = zeros(length(T),1);
sigma = 20/4.5; %currently, [sup(omega) in weeks]/4.5
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
eps = .1;
C0 = 2/eps; 
Cmax = C0*M.^(3/4); 
K = Cmax*eps/10;
mu0 = 0.05*M.^(-.25);
beta = 400/exp(1);
gma0 = ones(1,length(M))*5.*M.^(.75); %note minimum is about .12/eps
nu0 = ones(1,length(M))*5/(beta*eps).*M.^(-.25);

lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));
taubar = repmat(lambdab,length(T),1);

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000))));
lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));

R = exp(0:1:50)-1;
abc = [13, 23, 30];
ABC = R(abc);
%%
load('Nash0606.mat');

NashA = NASH0606{abc(1)};
[NA,~,~] = forwardTransport(T,M,omega,p,offset,R(abc(1)),Cmax,eps,K,mu0,gma0,nu0,NashA);
%
NashB = NASH0606{abc(2)};
[NB,~,~] = forwardTransport(T,M,omega,p,offset,R(abc(2)),Cmax,eps,K,mu0,gma0,nu0,NashB);
%
NashC = NASH0606{abc(3)};
[NC,~,~] = forwardTransport(T,M,omega,p,offset,R(abc(3)),Cmax,eps,K,mu0,gma0,nu0,NashC);
%%
figure(1)
SA = NA.*M.*M;
SB = NB.*M.*M;
SC = NC.*M.*M;
subplot(1,3,1)
semilogx(M,SA(end,:),'linewidth',1.2,'color',[.8,.8,0])
xlabel('m [mg]')
ylabel('S(50 weeks,m)')
title('Case A (\rho \approx 1.628\cdot 10^{5} 1/week)')
subplot(1,3,2)
semilogx(M,SB(end,:),'linewidth',1.2,'color',[0,.8,0])
xlabel('m [mg]')
axis([M(1) 1e6 0 70])
title('Case B (\rho = 3.585\cdot 10^{9} 1/week)')
subplot(1,3,3)
semilogx(M,SC(end,:),'linewidth',1.2,'color',[.8,0,.8])
xlabel('m [mg]')
title('Case C (\rho = 3.931\cdot 10^{12} 1/week)')
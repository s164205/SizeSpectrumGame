LSD = 1/17; %log-sate-difference
M = exp((1:200)*LSD)/exp(LSD);

dt = .03;
T = 1:dt:40;

offset = round(log(400)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
omega = zeros(length(T),1);
sigma = 10/4.5; %currently, [sup(omega) in weeks]/4.5
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
rho = 1e4;
eps = .2;
C0 = 1/eps; 
Cmax = C0*M.^(3/4); 
K = Cmax*eps/5;
mu0 = 0.01*M.^(-.25);
beta = 400/exp(1);
gma0 = ones(1,length(M))*13.*M.^(.75); %note minimum is about .12/eps
nu0 = ones(1,length(M))*13/(beta*eps).*M.^(-.25);

lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));
taubar = repmat(lambdab,length(T),1);

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000))));
lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));

R = exp(0:1:50)-1;
abc = [13, 20, 45];
ABC = R(abc);
%%
load('Real1306.mat');

taunash = Real1306{abc(1)};
[NA,gma,nu] = forwardTransport(T,M,omega,p,offset,R(abc(1)),Cmax,eps,K,mu0,gma0,nu0,taunash);
NashA = taunash;
%
taunash = Real1306{abc(2)};

[NB,gma,nu] = forwardTransport(T,M,omega,p,offset,R(abc(2)),Cmax,eps,K,mu0,gma0,nu0,taunash);
NashB = taunash;
%
taunash = Real1306{abc(3)};

[NC,gma,nu] = forwardTransport(T,M,omega,p,offset,R(abc(3)),Cmax,eps,K,mu0,gma0,nu0,taunash);
NashC = taunash;
%%
save('RCases.mat')
%%
figure(1)
subplot(1,3,1)
S = NA.*M.*M;

semilogx(M,S(end,:),'color',[0,0,.7])
axis([min(M) max(M) 0 max(S,[],'all')])
ylabel('S(m)')
xlabel('m [mg]')
legend('N(40 weeks,m)','location','sw')
title('\rho \approx 1.6\cdot 10^{5} 1/week')
%
subplot(1,3,2)
S = NB.*M.*M;

semilogx(M,S(end,:),'color',[0,0,.7])
axis([min(M) max(M) 0 max(S,[],'all')])
ylabel('S(m)')
xlabel('m [mg]')
legend('N(40 weeks,m)','location','sw')
title('\rho \approx 1.8\cdot 10^{8} 1/week')
%
subplot(1,3,3)
S = NC.*M.*M;

semilogx(M,S(end,:),'color',[0,0,.7])
axis([min(M) max(M) 0 max(S,[],'all')])
ylabel('S(m)')
xlabel('m [mg]')
legend('N(40 weeks,m)','location','sw')
title('\rho \approx 10^{19} 1/week')
axis([M(1) M(end) 0 .000025])
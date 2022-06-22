LSD = 1/17; %log-sate-difference
M = exp((1:200)*LSD)/exp(LSD);

dt = .03;
T = 1:dt:50;

offset = round(log(400)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
omega = zeros(length(T),1);
sigma = 20/4.5; 
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
rho = 1e5;
eps = .1;
C0 = 2/eps; 
Cmax = C0*M.^(3/4); 
K = Cmax*eps/10;
mu0 = 0.05*M.^(-.25);
beta = 400/exp(1);
gma0 = ones(1,length(M))*5.*M.^(.75);
nu0 = ones(1,length(M))*5/(beta*eps).*M.^(-.25);

lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));
taubar = repmat(lambdab,length(T),1);

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000))));
lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));
%%
[taunash,msqe,recordedTaubars,StepScalar] = searchBrouwers(T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar,20,1e-9,[],@bisectingH,8);

[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taunash);
[V,~] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
%%
save('InterestingNash.mat')
%%
figure(1)
set(1,'Position',[300 20 600 670])
subplot(2,2,1:2)
S = N.*M.*M;
J = [25,35,50];
P = [2.5,10,1474];
[~,Pint] = min(abs(P-M'));
for i = 1:numel(J)
    j = J(i);
    k = floor((j-1)/dt) +1;
    p = Pint(i);
    if floor(i/2) == i/2
       lty = ':'; 
    else
        lty = '-.';
    end
    semilogx(M,S(k,:),lty,'color',[0,0,.7])
    if i ~= 3
        txt = ['t = ' num2str(round(j))];
        text(M(p),S(k,p)+.01,txt,'VerticalAlignment','bottom')
    end
    axis([min(M) max(M) 0 max(S,[],'all')])
    hold on
end
semilogx(M,S(end,:),'-','color',[0,0,.7])
ylabel('S(m)')
text(M(p),S(k,p)+.003,['t = ' num2str(round(j)) ' -'],'HorizontalAlignment','right')
%
subplot(2,2,3:4)
semilogx(M,lambdab,'--','color',[.2,.2,.2])
text(4,lambdab(1)+.06,['\lambda_{b} \approx ' num2str(lambdab(1),2)])
hold on
J = [1,25,35,50];
P = [27,282,2222,1472];
[~,Pint] = min(abs(P-M'));
for i = 1:numel(J)
    j = J(i);
    k = floor((j-1)/dt) +1;
    p = Pint(i);
    if floor(i/2) == i/2
       lty = ':'; 
    else
        lty = '-.';
    end
    semilogx(M(1:end-1),taunash(k,1:end-1),lty,'color',[.7,0,0])
    
        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),taunash(k,p)+.01,txt,'HorizontalAlignment','right')
        
    axis([min(M) max(M) 0 1])
    hold on
end
semilogx(M(1:end-1),taunash(end,1:end-1),'-','color',[.7,0,0])
ylabel('\tau*(m)')
xlabel('m [mg]')
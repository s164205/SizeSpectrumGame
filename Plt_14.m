LSD = 1/17; %log-sate-difference
M = exp((1:200)*LSD)/exp(LSD);

dt = .01;
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
%%
[taunash,msqe,recordedTaubars,StepScalar] = searchBrouwers(T,M,VT,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar,20,1e-15,[],@bisectingH,8);
%taunash=taubar;
[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taunash);
[V,taunash] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
%[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taunash)
%%
save('NashRealistic.mat');
%%
figure(1)
set(1,'Position',[300 20 600 670])
subplot(2,2,1:2)
S = N.*M.*M;
J = [20,30,40];
P = [20,140,450];
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
    if i ~= 1
        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),S(k,p)+.003,txt,'HorizontalAlignment','right')
    end
    axis([min(M) max(M) 0 max(S,[],'all')])
    hold on
end
text(M(Pint(1)),S(floor((J(1)-1)/dt) +1,P(1))+ 3,['t = ' num2str(round(J(1)))],'VerticalAlignment','bottom')
semilogx(M,S(end,:),'-','color',[0,0,.7])
ylabel('S(m)')
%
subplot(2,2,3:4)
semilogx(M,lambdab,'--','color',[.2,.2,.2])
text(110,lambdab(1)-.05,['\lambda_{b} \approx ' num2str(lambdab(1),2)])
hold on
J = [1,20,30,40];
P = [110,45,390,1163];
[~,Pint] = min(abs(P-M'));
for i = 1:numel(J)
    j = J(i);
    k = floor((j-1)/dt) +1;
    p = Pint(i);
    if floor(i/3) == i/3
       lty = ':'; 
    elseif floor(i/3) == (i-1)/3
        lty = '--';
    else
        lty = '-.';
    end
    semilogx(M(1:end-1),taunash(k,1:end-1),lty,'color',[.7,0,0])
    
        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),taunash(k,p)+.01,txt,'HorizontalAlignment','right')
        
    axis([min(M) max(M) 0 1])
    hold on
end
%semilogx(M,taunash(105,:),'--','color',[.7,0,0])
%text(2.2,0.8,'-t = 1')
semilogx(M(1:end-1),taunash(end,1:end-1),'-','color',[.7,0,0])
ylabel('\tau*(m)')
xlabel('m [mg]')

%%
rho = 0;
[~,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
[V,taunash] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
[Nrho0,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taunash);
%%
Nsum = sum(N(end,:));
Ndiff = mean(abs(N(end,:)-Nrho0(end,:)),'all')./Nsum

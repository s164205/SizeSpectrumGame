tau = 0:.01:1;
gma = 5;
nu  = .5;
mu0 = .2;
cmax = 1;
K0 = .2;
g = @(tau) gma*tau./(1+gma*tau/cmax) - K0;
tauposg = tau(g(tau)>0);
mu = @(tau) nu*tau + mu0;
lambda = cmax*K0/(gma*(cmax-K0));

[~,maxdiff] = max(g(tau)-mu(tau));
%%
figure(1)
hold on
plot([0 1],[0 0],'k--')
plot([lambda lambda],[-K0,max(cmax,mu0+nu)],'--','color',[.7,.7,.7])
text(lambda,.6,['- \lambda = ' num2str(lambda)])
text(tau(maxdiff),.3,['- \tau* \approx ' num2str(tau(maxdiff))])
plot(tau,g(tau),'--','color',[0,.6,0],'linewidth',1.5)
p1 = plot(tauposg,g(tauposg),'color',[0,.6,0],'linewidth',1.5);
p2 = plot(tau,mu(tau),'color',[.9,.4,0],'linewidth',1.5);
plot([tau(maxdiff), tau(maxdiff)],[g(tau(maxdiff)),mu(tau(maxdiff))],'--','color',[.1,.1,.1])
plot([tau(maxdiff), tau(maxdiff)],[0,mu(tau(maxdiff))],'--','color',[.7,.7,.7])
axis([0 1 -K0 .8])
xlabel('\tau')
legend([p1 p2],{'$\frac{\partial V}{\partial m}\cdot g$','$\mu V$'},'location','se','interpreter','latex')
%%
LSD = 1/17; %log-sate-difference
M = exp((1:300)*LSD)/exp(LSD);
offset = round(log(500)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
taubar = ones(1,length(M))*1;
rho = 1e5;
nu0 = ones(1,length(M))*.15.*M.^(-.25);
gma0 = ones(1,length(M))*6.*M.^(.75);
N = zeros(1,length(M));
N(1,158) = 1e-3; %S(10^4) = 10

K = 0.2*M.^(.75);
mu0 = 0.05*M.^(-.25);

gma = growthVector(N,M,rho,gma0,taubar,p,offset);
nu = deathVector(N,rho,nu0,taubar,p,offset);
%%
figure(2)
loglog(M(1:end-1),gma(1:end-1)./M(1:end-1),'color',[0,.6,0],'linewidth',1.5)
hold on
loglog(M,K./M,'--','color',[0,.6,0],'linewidth',1.5)
loglog(M,nu,'color',[.9,.4,0],'linewidth',1.5)
loglog(M,mu0,'--','color',[.9,.4,0],'linewidth',1.5)
axis([M(1) M(end) 1e-4 1e2])
plot([M(157) M(157)],[1e-4 10],'color',[0,0,.7])
plot([M(157) M(158)],[10 10],'color',[0,0,.7])
plot([M(158) M(158)],[1e-4 10],'color',[0,0,.7])
legend({'$\gamma/m$','$(-K)/m$','$\nu$','$\mu_c$','$S$ [mg]'},'interpreter','latex')
ylabel('rate [1/week]')
xlabel('m [mg]')

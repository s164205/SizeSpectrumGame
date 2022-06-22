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
taubar = ones(length(T),length(M))*.8;
rho = 1e5;
eps = .1;
C0 = 2/eps; 
Cmax = C0*M.^(3/4); 
K = Cmax*eps/10;
mu0 = 0.05*M.^(-.25);
beta = 400/exp(1);
gma0 = ones(1,length(M))*5.*M.^(.75); %note minimum is about .12/eps
nu0 = ones(1,length(M))*5/(beta*eps).*M.^(-.25);

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000))));
lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));

[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
[V,tauopt] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
%%
figure(1)
set(1,'Position',[300 20 600 670])
subplot(2,2,1:2)
S = N.*M.*M;
J = [25,35,50];
P = [4,11,43,2650];
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

        txt = ['t = ' num2str(round(j)) ' - '];
        text(M(p),S(k,p)+.003,txt,'HorizontalAlignment','right')

    axis([min(M) max(M) 0 max(S,[],'all')])
    hold on
end
text(M(Pint(end))-1000,S(end,Pint(end)),'t = 50','VerticalAlignment','bottom')
semilogx(M,S(end,:),'-','color',[0,0,.7])
ylabel('S(m)')
%
subplot(2,2,3:4)
semilogx(M,lambdab,'--','color',[.2,.2,.2])
text(3,lambdab(1)+.06,['\lambda_{b} \approx ' num2str(lambdab(1),2)])
hold on
 J = [1,25,35,50];
P = [22,88,1100,1561];
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
    semilogx(M(1:end-1),tauopt(k,1:end-1),lty,'color',[.7,0,0])
    
        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),tauopt(k,p)+.01,txt,'HorizontalAlignment','right')
        
    axis([min(M) max(M) 0 1])
    hold on
end
semilogx(M(1:end-1),tauopt(end,1:end-1),'-','color',[.7,0,0])
ylabel('\tau*(m)')
xlabel('m [mg]')

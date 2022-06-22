LSD = 1/17; %log-sate-difference
M = exp((1:200)*LSD)/exp(LSD);

dt = .03;
T = 1:dt:50;

offset = round(log(500)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
omega = zeros(length(T),1);
sigma = 20/4.5; %currently, [sup(omega) in weeks]/4.5
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
taubar = ones(length(T),length(M))*1;
rho = 0;
eps = .6;
C0 = 2/eps; 
Cmax = C0*M.^(3/4); 
K = Cmax*eps/10;
mu0 = 0.05*M.^(-.25);
nu0 = ones(1,length(M))*1.*M.^(-.25);
gma0 = ones(1,length(M))*10/6.*M.^(.75);

lambda = Cmax.*K./(gma0.*(eps*Cmax-K));

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000))));

[~,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
[V,tauopt] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,tauopt);
%%
figure(1)
set(1,'Position',[300 20 600 670])
subplot(2,2,1:2)
J = [1,15,25,35,50];
P = [86000,86000,86000,86000,1400];
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
    semilogx(M,V(k,:),lty,'color',[.7,0,.7])

        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),V(k,p)+.01,txt,'HorizontalAlignment','right')

    axis([min(M) max(M) 0 1])
    hold on
end
semilogx(M,V(end,:),'-','color',[.7,0,.7])
ylabel('V(m)')
%
subplot(2,2,3:4)
semilogx(M,lambda,'--','color',[.2,.2,.2])
text(18,lambda(1)+.06,['\lambda \approx ' num2str(lambda(1),2)])
hold on
axis([min(M) max(M) 0 1])
P = [1.6,2.5,37.5,200,975];
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
    if i~=1 && i~=2
        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),tauopt(k,p)+.02,txt,'HorizontalAlignment','right')
    elseif i == 1
        txt = ['t = ' num2str(round(j))];
        text(M(p),tauopt(k,p)-.01,txt,'VerticalAlignment','top')
    else
        txt = ['t = ' num2str(round(j))];
        text(M(p),tauopt(k,p)-.01,txt,'VerticalAlignment','bottom')
    end
end
semilogx(M(1:end-1),tauopt(end,1:end-1),'-','color',[.7,0,0])
ylabel('\tau*(m)')
xlabel('m [mg]')
%%
S = N.*M.*M;
figure(2)
J = [15, 25,35,50];
P = [3,5,58,510];
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

        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),S(k,p),txt,'HorizontalAlignment','right')

    axis([min(M) max(M) 0 max(S,[],'all')])
    hold on
end
semilogx(M,S(end,:),'-','color',[0,0,.7])
ylabel('S(m)')
xlabel('m [mg]')
%

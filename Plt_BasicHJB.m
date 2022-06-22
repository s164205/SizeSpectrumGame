States = 1:200;
LSD = 1/17; %log-sate-difference
M = exp(States*LSD)/exp(LSD);

dt = .03;
T = 1:dt:50;

offset = 1;
p = 1;
omega = zeros(length(T),1);
sigma = 20/4.5; %currently, [sup(omega) in weeks]/4.5
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
taubar = ones(length(T),length(M));
rho = 0;
eps = 1;
C0 = 0; 
Cmax = ones(1,length(M))*Inf; 
K = 0;
mu0 = 0;
nu0 = ones(1,length(M))*.5.*M.^(-.25);
gma0 = ones(1,length(M))*.5.*M.^(.75);

VT = M>1e3;

[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
[V,tauopt] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
%%
figure(1)
set(1,'Position',[300 20 600 670])
subplot(4,4,1:12)
J = [1,30:5:50];
P = [5,125,250,427,685];
[~,Pint] = min(abs(P-M'));
for i = 1:numel(J)-1
    j = J(i);
    k = floor((j-1)/dt) +1;
    p = Pint(i);
    if floor(i/2) == i/2
       lty = ':'; 
    else
        lty = '-.';
    end
    semilogx(M,V(k,:),lty,'color',[.7,0,.7])
    if k ~= 1
        txt = ['t = ' num2str(round(j)) ' -'];
        text(M(p),V(k,p),txt,'HorizontalAlignment','right')
    end
    axis([min(M) max(M) 0 1])
    hold on
end
text(M(Pint(1)),V(1,Pint(1)),'t = 1','VerticalAlignment','bottom')
semilogx(M,V(end,:),'-','color',[.7,0,.7])
txt = ['t = ' num2str(round(T(end))) ' -'];
text(1e3,.75,txt,'HorizontalAlignment','right')
ylabel('V(m)')
%
subplot(4,4,13:16)
semilogx([1 1e3],[1, 1],[1e3 1e3],[0, 1],[1e3 1e6],[0, 0],'color',[.8,0,0],'linewidth',1.5)
ylabel('\tau^{*}(m)')
xlabel('m [mg]')
text(1e3,.5,'t \in [1,50] -','HorizontalAlignment','right')
axis([min(M) max(M) 0 1])
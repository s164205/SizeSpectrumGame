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



[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
%%
S = N.*M.*M;
figure(1)
J = 30:5:50;
P = [120,260,482,725,1160];
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
    semilogx(M,S(k,:),lty,'color',[0,0,.7])
    txt = ['-t = ' num2str(round(j))];
    text(M(p),S(k,p)+.3,txt,'HorizontalAlignment','left')
    axis([min(M) max(M) 0 1])
    hold on
end
i = numel(J);
p= Pint(i);
semilogx(M,S(end,:),'b-','linewidth',2)
xlabel('m [mg]')
ylabel('S(m)')
text(M(p),S(end,p)+.5,'-t = 50','HorizontalAlignment','left')
axis([min(M) max(M) 0 max(S,[],'all')])

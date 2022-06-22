LSD = 1/17; %log-sate-difference
M = exp((1:.5:200)*LSD)/exp(LSD);

dt = .03;
T = 1:dt:40;

offset = round(log(400)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
omega = zeros(length(T),1);
sigma = 10/4.5; %currently, [sup(omega) in weeks]/4.5
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
rho = 0;
eps = .2;
C0 = 1/eps; 
Cmax = C0*M.^(3/4); 
K = Cmax*eps/5;
mu0 = 0.01*M.^(-.25);
beta = 400/exp(1);
%gma0 = ones(1,length(M))*18.*M.^(.75); %note minimum is about .12/eps
%nu0 = ones(1,length(M))*18/(beta*eps).*M.^(-.25);

lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));
taubar = repmat(lambdab,length(T),1);

VT = ones(1,length(M))./(1+exp(-4*(log(M)-log(1000))));
lambdab = Cmax.*K./(gma0.*(eps*Cmax-K));
B = 1:40;
%%
taumean = zeros(length(B),1);
for i = 1:length(B)
    gma0 = ones(1,length(M))*B(i).*M.^(.75);
    nu0 = ones(1,length(M))*B(i)/(beta*eps).*M.^(-.25);

    [N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
    [V,taunash] = backwardHJB(T,M,VT,gma,nu,Cmax,eps,K,mu0);
    taumean(i) = mean(taunash,'all');
end
%%
figure(1)
plot(B,taumean,'.-','color',[.7,0,0])
hold on
plot([0 B(end)],[.5 .5],'--','color',[.2 .2 .2])
plot(13,taumean(13),'k.','markersize',7)
text(8.8,.48,'b = 13')
xlabel('b [1/week]')
ylabel('mean \tau*')
%% Constant biomass
M = 1:200;
dt = .01;
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
nu0 = ones(1,length(M))*5./M;
gma0 = ones(1,length(M))*5;


[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
%%
B = N.*M;
if ishandle(1)
    close 1
end
figure(1)
subplot(1,2,1)
hold on
J = 5:5:50;
for i = 1:numel(J)-1
    j = J(i);
    k = floor((j-1)/dt) +1;
    [modey,modex] = max(B(k,1:end-1));
    plot(B(k,1:end-1),'b-')
    text(modex,modey,['t = ' num2str(j)]);
end
[modey,modex] = max(B(end,1:end-1));
plot(B(end,1:end-1),'b-','linewidth',2)
text(modex,modey,'t = 50');
xlabel('m [mg]')
ylabel('B(m)')
subplot(1,2,2)
plot(T,sum(B,2)','linewidth',2)
xlabel('t')
ylabel('\int_{0}^{200}N(t,m)dm')
axis([1 T(end), 0 1.2])

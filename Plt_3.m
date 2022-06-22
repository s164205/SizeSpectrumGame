M = 1:200;
dt = .01;
T = 1:dt:50;

offset = 1;
p = 1;
omega = zeros(length(T),1);
sigma = 20/4.5; %currently, ([sup(omega) in weeks]/4.5)*sigma
Nomega = normpdf(-4:dt/sigma:4);
omega(1:length(Nomega)) = Nomega/sum(Nomega);
taubar = ones(length(T),length(M));
rho = 0;
eps = 1;
C0 = 0; 
Cmax = ones(1,length(M))*Inf; 
K = 0;
mu0 = 0;
nu0 = ones(1,length(M))*.02;
gma0 = ones(1,length(M))*3.5;


[N,gma,nu] = forwardTransport(T,M,omega,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
%% Number plot
figure(1)
hold on
J = 10:10:50;
P = [10,10,45,83,116];
for i = 1:numel(J)-1
    j = J(i);
    k = floor((j-1)/dt) +1;
    p = P(i);
    if floor(i/2) == i/2
       lty = ':'; 
    else
        lty = '-.';
    end
    plot(N(k,:),lty,'color',[0,0,.7])
    txt = ['-t = ' num2str(round(j))];
    text(p+.5,N(k,p)+.0003,txt,'HorizontalAlignment','left')
end
i = numel(J);
p= P(i);
plot(N(end,:),'b-','linewidth',2)
xlabel('m [mg]')
ylabel('N(m)')
text(p+.5,N(end,p)+.0003,'-t = 50','HorizontalAlignment','left')
%% Biomass plot
B = N.*M;
figure(1)
hold on
J = [15, 20:10:50];
P = [21,30,50,83,116];
for i = 1:numel(J)-1
    j = J(i);
    k = floor((j-1)/dt) +1;
    p = P(i);
    if floor(i/2) == i/2
       lty = ':'; 
    else
        lty = '-.';
    end
    plot(B(k,:),lty,'color',[0,0,.7])
    txt = ['-t = ' num2str(round(j))];
    text(p+.5,B(k,p)+.0003,txt,'HorizontalAlignment','left')
end
i = numel(J);
p= P(i);
plot(B(end,:),'b-','linewidth',2)
xlabel('m [mg]')
ylabel('B(m)')
text(p+.5,B(end,p)+.0003,'-t = 50','HorizontalAlignment','left')

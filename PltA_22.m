Dts = reshape((10.^(-(0:3))'*[5, 2 , 1])',[],1);

LSD = 1/17; %log-sate-difference
M = exp((1:200)*LSD)/exp(LSD);
offset = round(log(500)/LSD);
p = normpdf(-(offset-1)*LSD:LSD:(offset-1)*LSD);
N0 = [2*ones(floor(length(M)/2),1)/length(M); zeros(ceil(length(M)/2),1)];
taubar = ones(length(T),length(M))*1;
rho = 2e2;
eps = 1;
C0 = Inf; 
Cmax = Inf; 
K = 0;
mu0 = 0;
nu0 = zeros(length(M),1);
gma0 = nu0;
ActualMass = M*N0;

E = zeros(length(Dts),1);
for i = 1:length(Dts)
    disp(['starting test no. ' num2str(i) ' with timesteps ' num2str(Dts(i))])
    
    T = 0:Dts(i):100;
    taubar = ones(length(T),length(M));
    N = forwardTransportWithN0(T,M,N0,p,offset,rho,Cmax,eps,K,mu0,gma0,nu0,taubar);
    totMass = M*N(end,:)';
    E(i) = abs(totMass-ActualMass);
end
save('E.mat','E')
%%
figure(1)
subplot(2,1,1)
plot(E,'.-','linewidth',1.5,'markersize',10)
ylabel('E^{2}')
axis([0 length(E) 0 5])
strs = ['\Detla t = ', num2str(Dts(end))];
for i = 1:length(Dts)-1
    strs = char(['\Delta t = ', num2str(Dts(end-i))],strs);
end
xticklabels(strs)
subplot(2,1,2)
loglog(Dts,E,'.-','linewidth',1.5,'markersize',10)
grid on
ylabel('E^{2}')
xlabel('\Delta t')

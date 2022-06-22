load('InterestingNash.mat')
tauInt = taunash;
NInt = N;
TInt = T;
lambInt = lambdab;
load('NashRealistic.mat')
tauReal = taunash;
NReal = N;
TReal = T;
lambReal = lambdab;
%%
meanI = mean(tauInt,2);
meanR = mean(tauReal,2);
%%
figure(1)
plot(TReal,meanR,'linewidth',1.5)
hold on
plot(TInt,meanI,'linewidth',1.5)
legend({'Probable Parameters','Feasible Parameters'},'location','se')
xlabel('t [weeks]')
ylabel('mean \tau*(t)')
axis([1 50 0 1])

load('RealisticNash.mat');
Vr = V;
Tr = T;
load('InterestingNash.mat');
Vi = V;
Ti = T;
%%
Fomr = Vr(:,3); %odd things happen before m = 3
Fomi = Vi(:,3);
%%
figure(1)
%subplot(1,2,1)
plot(Tr,Fomr,'linewidth',1.5)
hold on
plot(Ti,Fomi,'linewidth',1.5)
axis([1 20 0 5e-3])
xlabel('t [weeks]')
ylabel('V(t,1mg)')
legend('Probable Parameters','Feasible Parameters')
% subplot(1,2,2)
% semilogy(Tr,Fomr)
% hold on
% semilogy(Ti,Fomi)
% axis([1 50 1e-15 1])
% xlabel('t [weeks]')
% ylabel('F^{(Spawn Time)}(t)')
% legend('Probable Parameters','Realistic Parameters')

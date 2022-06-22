load('NashInteresting')
FInt = N(end,:)*VT';
load('FPhi.mat')
load('INFOPHI.mat') %loading data from S_PHI9 and S_B
load('FB.mat')
load('INFOB.mat')
Phi0 = 1./PHI0; %I originally did it wrong
%%
DefPhi0 = 1/10;
DefB = 5;
figure(1)
subplot(1,2,1)
semilogx([1,Phi0,1e-5],[0,FPhi'],'linewidth',1.5)
xlabel('\phi_{0}')
ylabel('F')
axis([1e-5 1 0 1e-4])
subplot(1,2,2)
semilogx([1,B],[0,FB'],'linewidth',1.5)
xlabel('b [1/week]')
axis([1 max(B) 0 1e-4])
load('F0606.mat');
Fi = F0606;
load('FR1306.mat');
Fr = FR1306;
load('PARMSR1306.mat'); %First load data from S_recruitment
%%
R = exp(0:1:50)-1;
abc = [13, 23, 30];
ABC = R(abc);
%%
figure(1)
loglog(R,Fr,'.-',R,Fi,'.-','linewidth',1.5,'markersize',8)
hold on
plot(R(11),Fi(11),'o','markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[.8,.8,.8])
plot(ABC(1),Fi(abc(1)),'o','markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[.8,.8,0])
plot(ABC(2),Fi(abc(2)),'o','markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0,.8,0])
plot(ABC(3),Fi(abc(3)),'o','markersize',5,'MarkerEdgeColor','k','MarkerFaceColor',[.8,0,.8])
xlabel('\rho [1/week]')
ylabel('F')
legend({'Probable Parameters','Feasible Parameters','\rho \approx 10^{5} 1/week','Case A (\rho \approx 1.5\cdot 10^{5})',...
    'Case B (\rho \approx 3.6\cdot 10^{9})','Case C (\rho \approx 4\cdot 10^{12})'},'location','sw')
axis([.5 R(end)*2 1e-15 1])
%%
figure(2)
p3 = loglog(R,Fr(1).*R','--','color',[0.1 0.2470 0.5410]);
hold on
p4 = loglog(R,Fi(1).*R','--','color',[0.6500 0.1250 0.1980]);
p1 = loglog(R,Fr.*R','.-','linewidth',1.5,'markersize',8,'color',[0 0.4470 0.7410]);
p2 = loglog(R,Fi.*R','.-','linewidth',1.5,'markersize',8,'color',[0.8500 0.3250 0.0980]);
xlabel('$\overline{\Omega}$','interpreter','latex')
ylabel('$F\cdot \overline{\Omega}$','interpreter','latex')
legend([p1 p2 p3 p4],{'Realistic Parameters','Feasible Parameters',...
    '$F(\rho = 0)\cdot \overline{\Omega} (Real. parm.)$','$F(\rho = 0)\cdot \overline{\Omega} (Feas. parm.)$'}...
    ,'interpreter','latex','location','nw')
axis([1 R(end)*1.1 1e-10 1e20])

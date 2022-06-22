tauR = 0:.01:1;
tauI = 0:.01:1;
[TAUR,TAUI] = meshgrid(tauR,tauI);

F = @(tau) -(tau-.5).^2;

diffFtau = @(tauR,tauI) F(tauI)-F(tauR);

tr = [.3,.3];
tia = [.3,.2];
tib = [.3,.4];
trfut = [.4,.4];


%%
figure(1)
hold on
patch([0 0 .5],[0 1 .5],[0,.7,0])
patch([1 1 .5],[0 1 .5],[0,.7,0])
plot(tauR,tauI,'k','linewidth',2,'color',[0,0,.5])
contour(tauR,tauI,diffFtau(TAUR,TAUI),-.3:.1:.3,'showtext','on','color','black')
xlabel('\tau_{R}')
ylabel('\tau_{I}')
text(.2,.5,'+','fontsize',25)
text(.75,.5,'+','fontsize',25)
text(.485,.23,'-','fontsize',25)
text(.485,.8,'-','fontsize',25)
plot(tr(1),tr(2),'.k',tia(1),tia(2),'.r',tib(1),tib(2),'.k',trfut(1),trfut(2),'.k','markersize',20)
text(tr(1)+.02,tr(2),'\tau_{R,now}')
text(tia(1)+.02,tia(2),'\tau_{I,a}')
text(tib(1)+.01,tib(2)+.04,'\tau_{I,b}')
text(trfut(1)+.02,trfut(2),'\tau_{R,future}')
plot([tr(1),tia(1)],[tr(2),tia(2)],'k--')
plot([tr(1),tib(1)],[tr(2),tib(2)],'k--')
plot([trfut(1),tib(1)],[trfut(2),tib(2)],'k--')
arrow = annotation('textarrow',[tr(1),(tr(1)+trfut(1))*.5]+.015,[tr(2),(tr(2)+trfut(2))*.5]+.015,'FontSize',3,'Linewidth',2,'color',[0,0,.5]);
arrow.Parent = gca;
patch([0 0 .9 .9],[0 0 .9 .9],'white')
%annotation('textbox', [0.8, 0.9, 0.9, 0.9], 'String','$F(\tau_I)-F(\tau_R)$','interpreter','latex');
%annotation('textbox', [0.8, 0.8, 0.9, 0.9], 'String','$F(\tau) = (\frac{1}{2}-\tau)^2$','interpreter','latex');

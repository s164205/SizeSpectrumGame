xx = 0:.01:1; %resident
yy = 0:.01:1; %invader

[X,Y] = meshgrid(xx,yy);

f = @(x,y) y.*(1-2*x.^2);

F = f(X,Y);
xxopt = 1/sqrt(6);
yyopt = 1/sqrt(2);


%%
figure(1)
contour(xx,yy,F,'showtext','on')
hold on
p1 = plot(xx,yy,'linewidth',2);
p2 = plot(xx(xx<=yyopt),ones(sum(xx<=yyopt),1),'m','linewidth',2);
plot(xx(xx>=yyopt),zeros(sum(xx>=yyopt),1),'m','linewidth',2)
plot([yyopt yyopt],[0 1],'m','linewidth',2)
p3 = plot(xxopt,xxopt,'k.','markersize',20);
p4 = plot(yyopt,yyopt,'r.','markersize',20);
axis([0 1 0 1])
legend([p1,p2,p3,p4],{'$\tau = \overline{\tau}$','$\tau^{*}(\overline{\tau})$','$\tau^* (mono)$','$\tau^* (poly)$'},'interpreter','latex')
xlabel('$\overline{\tau}$','interpreter','latex')
ylabel('$\tau$','interpreter','latex')

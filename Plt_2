M = 1.25:.5:6.75;
Mintegrated = M-.25;
x = normpdf(M-4);
x = x/sum(x);
delta = M(2)-M(1);
bounds = [(M-delta/2) M(end)+delta/2];

Mprecise = .9:.01:7.1;
xprecise = normpdf(Mprecise-4);

%Note that this plot is TECHNICALLY incorrect, as the histogram doesn't
%actually show the integral of the precise distribution, but rather its
%value at the mid-point.

figure('DefaultAxesFontSize',14)
set(gcf, 'Position', get(0, 'Screensize'));
subplot(1,2,1)
bar(M,x)
hold on
for i = 1:length(bounds)
    plot([bounds(i), bounds(i)],[0 .5],'--','color',[.5 .5 .5])
end
axis([.9 7.1 0 .5]);
xlabel('m [mg]')
ylabel('P(m_{i}<N<m_{i+1})')
lgnd1 = legend('Fraction of population between m_{i} and m_{i+1}','indices m_{i}');
set(lgnd1,'FontSize',15);
subplot(1,2,2)
p2 = plot(Mprecise,xprecise,'k--','linewidth',.5);
axis([.9 7.1 0 .5]);
hold on
p1 = plot(Mintegrated,x/delta,'.-','linewidth',2,'markersize',17);
xlabel('m [mg]')
ylabel('pdf(N)')
lgnd2 = legend([p1,p2],{'Number distribution','Theoretical probability density function'});
set(lgnd2,'FontSize',15);


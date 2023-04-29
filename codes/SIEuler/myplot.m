function res = myplot(time,xb,x,rhox,ux,pb,rhoE,auxpar)


% if (plotcall == 0)
%     fig1 = figure( 'Position', [2 2 927 994]);
%     plotcall = 1;
% end

fontsize = 14;

title(strcat('time = ',num2str(time,3),', residual =',num2str(auxpar.res),'Interpreter','latex'))

subplot(1,3,1)
% plot(x,rho,'o-','MarkerFaceColor',[1 0 0])
plot(x,rhox,'-','Color',auxpar.linecolor)
grid on
set(gca,'GridLineStyle',':')
title({'$\rho$'},'Interpreter','latex','FontSize',fontsize)


subplot(1,3,2)
% plot(xv,vv,'o-','MarkerFaceColor',[0 1 0])
plot(x,ux,'-','Color',auxpar.linecolor)
grid on
set(gca,'GridLineStyle',':')
title({'$u$'},'Interpreter','latex','FontSize',fontsize)

subplot(1,3,3)
% plot(x,a,'o-','MarkerFaceColor',[0 0.5 0.5])
plot(xb,pb,'-','Color',auxpar.linecolor)
grid on
set(gca,'GridLineStyle',':')
title('$p$','Interpreter','latex','FontSize',fontsize)
pause(0.01)

%%
cd /Users/npizzo/Documents/Research/Wave_vortex/
ao=(1/2)^(2/3);
t=0:.01:2*pi;
clf
p1=plot(-ao*sin(t),ao*cos(t),'k','linewidth',2);
hold on
plot(-ao*sin(0),ao*cos(0),'ok','markersize',8,...
       'markerfacecolor','k')
hold on
p2=plot(cos(t)-ao*sin(t),ao*cos(t)+sin(t),'r', 'linewidth',2);
hold on
plot(cos(0)-ao*sin(0),ao*cos(0)-sin(0),'or',...
       'markerfacecolor','r','markersize',8)
set(gca,'fontsize',22)
xlabel('x position','interpreter','latex')
ylabel('y position','interpreter','latex')
l1=legend([ p1 p2],'Vortex','Wave packet')
set(l1,'interpreter','latex')

%% 
Lo=0.5*(-3+2*sqrt(2)); % stable
% Lo=-2.9651;
% Lo=0.5*(5+2*sqrt(2)); %unstable 
% Lo=0.5*(-1-2*sqrt(2));
Ho=(1+2*Lo)/(-1+2*Lo); 
a=sqrt((Lo-1/2*1^2)/(Ho-sqrt(1)));
k=0:0.01:20;
E1=k.^2.*(sqrt(k) - Ho).^2./((Lo - 1/2*k.^2).^2);
E2=(1-2./k.^2.*(Lo-0.5*k.^2).*(Ho-sqrt(k))).^2;
E3=E1.*(-1+E2);
% gE3=gradient(E3)./gradient(k);
% ggE3=gradient(gE3)./gradient(k);
clf
% subplot(2,1,1)
plot(k,E1.*(-1+E2),'k','linewidth',2)
hold on
% plot(k,ggE3,'r')
ylim([-1 10])
xlim([0 7])
set(gca,'fontsize',22)
xlabel('$\rho$','interpreter','latex')
ylabel('$\Pi(\rho)$','interpreter','latex')
set(gca,'YTickLabel',[])
% l1=legend('$\Pi(\rho)$','$\Pi''(\rho)$');
% set(l1,'interpreter','latex')
% title('$\mathcal{L}_0=1/2(-3-2\sqrt{2})$','interpreter','latex')
grid on
%% solve system numerically 
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,y] = ode45(@W1V1,[0:.05:12],[1+0.01;0.0;-0.01;a],options);
% check conservation laws
H=(y(:,1).^2+y(:,2).^2).^(1/4)-(y(:,1).*y(:,4)-y(:,2).*y(:,3))./...
    (y(:,3).^2+y(:,4).^2);
L=1/2*(y(:,1).^2+y(:,2).^2)-(y(:,1).*y(:,4)-y(:,2).*y(:,3));
% theta=atan(y(:,2)./y(:,1));
% v=y;
clf
l1=plot(x_sp,y_sp,'b','linewidth',1/2);
hold on
l2=plot(x_sv,y_sv,'c','linewidth',1/2);
hold on
l3=plot(y(:,2),-y(:,1),'--r','linewidth',1/8);
hold on
l4=plot(y(:,3)+y(:,2),y(:,4)-y(:,1),'--k','linewidth',1/8);
hold off
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
r=1.5;
xlim([-r r])
ylim([-r r])
l3=legend('Vortex','Wave packet');
set(l3,'interpreter','latex')
% figure
% subplot(2,1,1)
% plot(t,y(:,1),'--r');
% hold on
% plot(t,y(:,2),'--k');
% hold on
% plot(t,y(:,3),'--b');
% hold on
% plot(t,y(:,4),'--c');
% subplot(2,1,2)
% plot(t,L)
% hold on
% plot(t,H,'k')
% % compute rho from the numerics
% k_n=sqrt(y(:,1).^2+y(:,2).^2);
% E1_n=k_n.^2.*(sqrt(k_n) - H(1)).^2./((L(1) - 1/2*k_n.^2).^2);
% E2_n=(1-2./k_n.^2.*(L(1)-0.5*k_n.^2).*(H(1)-sqrt(k_n))).^2;
% E3_n=E1_n.*(-1+E2_n);
% k=0:0.01:10;
% E1=k.^2.*(sqrt(k) - H).^2./((L - 1/2*k.^2).^2);
% E2=(1-2./k.^2.*(L-0.5*k.^2).*(H-sqrt(k))).^2;
% E3=E1.*(-1+E2);
% % figure 
% plot(k,E3)
% ylim([-10 10])

%% animate this? 
clf
for i=1:length(t)
%     subplot(2,1,1) 
plot(y(i,2),-y(i,1),'or',...
    'markerfacecolor','r','markersize',8);
hold on
plot(y(:,2),-y(:,1),'--r','linewidth',1/8);
hold on
plot(y(i,3)+y(i,2),y(i,4)-y(i,1),'ok',...
    'markerfacecolor','k','markersize',8);
hold on
plot(y(:,3)+y(:,2),y(:,4)-y(:,1),'--k','linewidth',1/8);
hold off
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
r=3;
xlim([-r r])
ylim([-r r])
% subplot(2,1,2)
% % plot(k,E1.*(-1+E2),'k','linewidth',2)
% plot(k,E3,'k','linewidth',2)
% hold on
% plot(sqrt(y(i,1).^2+y(i,2).^2),E3_n(i,1),'ro',...
%     'markerfacecolor','r','markersize',8)
% hold off
% xlim([0 6])
% ylim([-1 3])
% set(gca,'fontsize',22)
% xlabel('$\rho$','interpreter','latex')
% ylabel('$\Pi(\rho)$','interpreter','latex')
pause(0.01)
end
%% save as a movie
% video settings
frames = 2000;
vfrate = 10; % frame rate of video
vname = 'bound_2';
fhand = figure('visible','off');
set(gcf,'color','w');
set(gcf,'Renderer','zbuffer');
writerObj = VideoWriter(vname,'MPEG-4');
writerObj.FrameRate = vfrate;
writerObj.Quality = 100;
open(writerObj);
clf
for i=1:3:length(t)
    subplot(2,1,1)
plot(y(i,2),-y(i,1),'or',...
    'markerfacecolor','r','markersize',8);
hold on
plot(y(:,2),-y(:,1),'--r','linewidth',1/8);
hold on
plot(y(i,3)+y(i,2),y(i,4)-y(i,1),'ok',...
    'markerfacecolor','k','markersize',8);
hold on
plot(y(:,3)+y(:,2),y(:,4)-y(:,1),'--k','linewidth',1/8);
hold off
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
r=3;
xlim([-r r])
ylim([-r r])
pbaspect([1 1 1])
subplot(2,1,2)
% plot(k,E1.*(-1+E2),'k','linewidth',2)
plot(k,E3,'k','linewidth',2)
hold on
plot(sqrt(y(i,1).^2+y(i,2).^2),E3_n(i,1),'ro',...
    'markerfacecolor','r','markersize',8)
hold off
xlim([0 6])
ylim([-1 3])
set(gca,'fontsize',22)
xlabel('$\rho$','interpreter','latex')
ylabel('$\Pi(\rho)$','interpreter','latex')
    drawnow;
    cframe = getframe(gcf); % take shot
    writeVideo(writerObj,cframe); % add it to video
end
close(writerObj); 
%% plot orbits
clf
l1=plot(x_sp,y_sp,'b','linewidth',1/2);
hold on
l2=plot(x_sv,y_sv,'c','linewidth',1/2);
hold on
l3=plot(y(:,2),-y(:,1),'--r','linewidth',1/8);
hold on
l4=plot(y(:,3)+y(:,2),y(:,4)-y(:,1),'--k','linewidth',1/8);
hold off
set(gca,'fontsize',12)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
r=1.5;
xlim([-r r])
ylim([-r r])
l3=legend('Vortex','Wave packet','Vortex (Perturbed)',...
    'Wave packet (Perturbed)');
set(l3,'interpreter','latex')

%% roots of unity plot 
N=15;
clf
l1=plot(0,0,'or','markerfacecolor','r');
hold on
l2=plot(real(exp(2*pi*I*(1:N)/N)),...
    imag(exp(2*pi*I*(1:N)/N)),'ok','markerfacecolor','k');
hold off
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
r=1.5;
xlim([-r r])
ylim([-r r])
l3=legend('Vortex','Wave packet');
set(l3,'interpreter','latex')
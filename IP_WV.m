% Scratch for solving for N vortices and M wave packets 
cd /Users/npizzo/Documents/Research/Wave_vortex
tstart=cputime;
N = 15;
M = 1;
L = 3;
init=zeros(4*M+2*N,1);
init(M+1:2*M,1)=2*pi/M*[1:M];
init(2*M+1:3*M,1)=ones(length(M),1);
% init(3*M+1:4*M,1)=0.1*(rand(M,1)-1/2);
init(3*M+1:4*M,1)=zeros(M,1);
init(4*M+1:4*M+N,1)=2*pi*rand(N,1);
% init(4*M+N+1:4*M+2*N,1)=1/2+[1:N];
init(4*M+N+1:4*M+2*N,1)=2*pi*rand(N,1);
A=1*ones(M,1); 
G=0.01*(ones(N,1));
% G=s( 2*rem(1:N,2) - 1);
G(1:3,1)=0.01*2*pi*ones(3,1);
g = 1; 
options = odeset('RelTol',1e-4,'AbsTol',1e-5);
[T,Y] = ode15s(@(t,y)WMVN_IP( t, y, N, M, G, A, g, L), [0:0.1:2000], init, options);
tend=cputime-tstart
% Y(15,2)
%% display data
clf 
plot(mod(Y(:,1:M),2*pi), mod(Y(:,M+[1:M]),2*pi), 'bd',...
    'markerfacecolor','b');
hold on
p0=plot(mod(Y(:,1),2*pi), mod(Y(:,M+1),2*pi), 'bd',...
    'markerfacecolor','b');
hold on
plot(mod(Y(:,4*M+[4:N]),2*pi), mod(Y(:,4*M+N+[4:N]),2*pi),'r.',...
    'markerfacecolor','r');
p1=plot(mod(Y(:,4*M+4),2*pi), mod(Y(:,4*M+N+4),2*pi),'r.',...
    'markerfacecolor','r');
hold on
plot(mod(Y(:,4*M+[1:3]),2*pi), mod(Y(:,4*M+N+[1:3]),2*pi),'ko',...
    'markerfacecolor','k');
p2=plot(mod(Y(:,4*M+1),2*pi), mod(Y(:,4*M+N+1),2*pi),'ko',...
    'markerfacecolor','k');
hold on
% plot(mod(Y(end,1:M),2*pi), mod(Y(end,M+[1:M]),2*pi), 'co');
% hold on 
% plot(mod(Y(end,4*M+[1:N]),2*pi), mod(Y(end,4*M+N+[1:N]),2*pi),'bo',...
%     'markerfacecolor','b');
xlim([0 2*pi])
ylim([0 2*pi]) 
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
legend([p1 p2 p0],...
    'Weak vortices', 'Strong vortices','Wave packet')
% pause(0.05)
%% histogram this is not the best way to show clumpiness
% make one long vector
Vx=zeros(length(T)*N,1);
Vy=Vx;
Lt=length(T);
for i=1:N
    Vx(1+(i-1)*Lt:Lt*(1+(i-1)),1)=...
        mod(Y(:,4*M+i),2*pi);
    Vy(1+(i-1)*Lt:Lt*(1+(i-1)),1)=...
        mod(Y(:,4*M+N+i),2*pi);
end
%
clf
nbins= [10 10];
[NN,C]=hist3([Vx, Vy],nbins);
contourf(C{1},C{2},NN')
colorbar
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
title('Histogram: 10 vortices, 1 wave packet','interpreter','latex')
%%
clf
plot(T,sqrt(Y(:,2*M+1).^2+Y(:,3*M+1).^2),'r')
set(gca,'fontsize',22)
xlabel('$t$','interpreter','latex')
ylabel('$|\bf{k}|$','interpreter','latex')
%% movie
clf
for i=1:10:length(T)
plot(mod(Y(i,1:M),2*pi),mod(Y(i,M+[1:M]),2*pi),'*','Color',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)],'markerfacecolor',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)]);
hold on
plot(mod(Y(i,4*M+[1:N]),2*pi),mod(Y(i,4*M+N+[1:N]),2*pi),'o','Color',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)]);
xlim([0 2*pi])
ylim([0 2*pi])
set(gca,'fontsize',22)
xlabel('X','interpreter','latex')
ylabel('Y','interpreter','latex')
pause(0.1)
end
%%
clf 
for j=1:length(T)
p1=plot(Y(j,1),Y(j,M+1),'o','Color',...
    [1-j/numel(T),1-j/numel(T),1-j/numel(T)]);
hold on
p2=plot(Y(j,2),Y(j,M+2),'d','Color',...
    [1-j/numel(T),1-j/numel(T),1-j/numel(T)]);
hold on
p3=plot(Y(j,4*M+1),Y(j,4*M+N+1),'*','Color',...
    [1-j/numel(T),1-j/numel(T),1-j/numel(T)]);
hold on
p4=plot(Y(j,4*M+2),Y(j,4*M+N+2),'*','Color',...
    [1-j/numel(T),1-j/numel(T),1-j/numel(T)]);
hold on
p3=plot(Y(j,4*M+3),Y(j,4*M+N+3),'*','Color',...
    [1-j/numel(T),1-j/numel(T),1-j/numel(T)]);
hold on
p4=plot(Y(j,4*M+4),Y(j,4*M+N+4),'*','Color',...
    [1-j/numel(T),1-j/numel(T),1-j/numel(T)]);
hold on
xlim([-10 10])
ylim([-10 10])
pause(0.05)
end
% set(gca,'fontsize',22)
% xlabel('X','interpreter','latex')
% ylabel('Y','interpreter','latex')
% l1=legend([p1 p2 p3],'Wave packet 1','Wave packet 2','Vortex');
% set(l1,'interpreter','latex');
% colormap(flipud(gray(30)));
% cbh=colorbar;
% set(cbh,'YTick',[0:.1:T(30)])
% set(get(cbh,'label'),'string','time')
%%
% check hamiltonian conservation 
H=zeros(length(T),1);
for j=1:length(T)  
    s=0;     S1=0;    r=0;
    for k=1:length(G)  
    for p=1:length(A)
        for n=-L:L
            for m=-L:L
     s=s+...
    (A(p)*G(k)/2/pi*...
    ((Y(j,p)-Y(j,4*M+k)+2*pi*n)*Y(j,3*M+p)-...
    (Y(j,M+p)-Y(j,4*M+N+k)+2*pi*m)*Y(j,2*M+p)))/...
    ((Y(j,p)-Y(j,4*M+k)+2*pi*n).^2+...
    (Y(j,p+M)-Y(j,4*M+N+k)+2*pi*m).^2);
            end
        end
    end
    end
    for p=1:length(A)
    S1=S1+(A(p).*sqrt(g*...
    sqrt(Y(j,2*M+p).^2+Y(j,3*M+p).^2)));
    end
    for k=1:length(G)
    if numel(G)>1
    for l=k+1:length(G)
        for n=-L:L
            for m=-L:L
      r=r+G(k)*G(l)*log(((Y(j,4*M+k)-Y(j,4*M+l)+2*pi*n).^2+...
          (Y(j,4*M+N+k)-Y(j,4*M+N+l)+2*pi*m).^2));
            end
        end
    end
    else
    r=0;
    end
    end
    H(j,1)=S1-1/2/pi*r+s;
end
clf
plot(T,H/H(1))
set(gca,'fontsize',20)
xlabel('t','interpreter','latex')
ylabel('H','interpreter','latex')
% ylim([0 2])
%% momentum
P=zeros(length(T),2);
for j=1:length(T)  
    S1=[0 0];    r=[0 0];
    for p=1:length(A)
    S1=S1+A(p).*...
    [Y(j,2*M+p),Y(j,3*M+p)];
    end
    for k=1:length(G)
        for n=-L:L
            for m=-L:L
      r=r+G(k)*[Y(j,4*M+N+k)+2*pi*m,-Y(j,4*M+k)+2*pi*n];
            end
        end
    end
    P(j,:)=S1+r;
end
%
clf
plot(T,P(:,1)/P(1,1),'r')
hold on
plot(T,P(:,2)/P(1,2),'b')
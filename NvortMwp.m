% Scratch for solving for N vortices and M wave packets 
cd /Users/npizzo/Documents/Research/Wave_vortex
N = 1;
M = 5;
% Case 1: 2 wave packets, 1 vortex, circular motion
% r = 1/16;
% chi = 1;
% phi = 0;
% init = [chi*cos(phi); -chi*cos(phi); ...
%     chi*sin(phi); -chi*sin(phi); r*sin(phi); -r*sin(phi);...
%    -r*cos(phi); r*cos(phi); 0; 0;];
% A = - [1 1];
% G = 2 * pi; 
% 
init=zeros(4*M+2*N,1);
init(M+1:2*M,1)=[1:M];
init(2*M+1:3*M,1)=ones(length(M),1);
% init(3*M+1:4*M,1)=0.1*(rand(M,1)-1/2);
init(3*M+1:4*M,1)=zeros(M,1);
init(4*M+1:4*M+N,1)=2+0*[1:N]'+zeros(N,1);
% init(4*M+N+1:4*M+2*N,1)=1/2+[1:N];
init(4*M+N+1:4*M+2*N,1)=1/4+0.05*[1:N];
A=.1*ones(M,1); 
% G=.01*([1:N]);
G=( 2*rem(1:N,2) - 1);
% init=10*rand(4*M+2*N,1)-5;
% init=zeros(4*M+2*N,1);
% init(M+1:2*M,1)=[1:M];
% init(2*M+1:3*M,1)=ones(length(M),1);
% init(3*M+1:4*M,1)=0.1*(rand(M,1)-1/2);
% init(3*M+1:4*M,1)=zeros(M,1);
% init(4*M+1:4*M+N,1)=3+0*[1:N]'+zeros(N,1);
% init(4*M+N+1:4*M+2*N,1)=1/2+[1:N];
% init(4*M+N+1:4*M+2*N,1)=4+1/4+0.05*[1:N];
% A=.1*ones(M,1); 

% G=.01*( 2*rem(1:N,2) - 1);
% G=.01*([1:N]);
g = 1;
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[T,Y] = ode45(@(t,y)WMVN( t, y, N, M, G, A, g), [0:.1:100], init, options);
% Y(15,2)
%% display data

clf 
plot( Y(:,1:M), Y(:,M+[1:M]), 'ko');
hold on
plot( Y(:,4*M+[1:N]), Y(:,4*M+N+[1:N]),'rd');
% xlim([-10 10])
% ylim([-10 10])
% pause(0.05)
%% movie
clf
for i=1:length(T)
plot(Y(i,1:M),Y(i,M+[1:M]),'o','Color',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)],'markerfacecolor',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)]);
hold on
plot(Y(i,4*M+[1:N]),Y(i,4*M+N+[1:N]),'*','Color',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)]);
xlim([-1 10])
ylim([0 20])
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
    s=0;     S1=0;    L=0;
    for k=1:length(G)  
    for p=1:length(A)
     s=s+...
    (A(p)*G(k)/2/pi*...
    ((Y(j,p)-Y(j,4*M+k))*Y(j,3*M+p)-...
    (Y(j,M+p)-Y(j,4*M+N+k))*Y(j,2*M+p)))/...
    ((Y(j,p)-Y(j,4*M+k)).^2+...
    (Y(j,p+M)-Y(j,4*M+N+k)).^2);
    end
    end
    for p=1:length(A)
    S1=S1+(A(p).*sqrt(g*...
    sqrt(Y(j,2*M+p).^2+Y(j,3*M+p).^2)));
    end
    for k=1:length(G)
    if numel(G)>1
    for l=k+1:length(G)
      L=L+G(k)*G(l)*log(((Y(j,4*M+k)-Y(j,4*M+l)).^2+...
          (Y(j,4*M+N+k)-Y(j,4*M+N+l)).^2));
    end
    else
    L=0;
    end
    end
    H(j,1)=S1-1/2/pi*L+s;
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
      r=r+G(k)*[Y(j,4*M+N+k),-Y(j,4*M+k)];
    end
    P(j,:)=S1+r;
end
%
clf
subplot(2,1,1)
plot(T,P(:,1),'r')
subplot(2,1,2)
plot(T,P(:,2),'b')
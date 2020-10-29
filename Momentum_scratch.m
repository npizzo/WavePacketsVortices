%% correcting the momentum
% N = 15;
cd /Users/npizzo/Documents/Research/Wave_vortex
tstart=cputime;
N = 10; 
M = 0;
L = 3;
init=zeros(4*M+2*N,1);
init(1:M,1)=2*pi*rand(M,1);
init(M+1:2*M,1)=2*pi*rand(M,1);
init(2*M+1:3*M,1)=ones(length(M),1);
init(3*M+1:4*M,1)=zeros(M,1);
init(4*M+1:4*M+N,1)=rand(N,1);
init(4*M+N+1:4*M+2*N,1)=rand(N,1);
A=1*ones(M,1); 
G=(ones(N,1));
g = 1; 
%
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[T,Y] = ode15s(@(t,y)WMVN_IP( t, y, N, M, G, A, g, L), [0:.1:50], init, options);
tend=cputime-tstart
% Y(15,2)
%% display data
% clf 
figure
plot(mod(Y(:,1:M),2*pi), mod(Y(:,M+[1:M]),2*pi), 'bd',...
    'markerfacecolor','b');
hold on
plot(mod(Y(:,4*M+[1:N]),2*pi), mod(Y(:,4*M+N+[1:N]),2*pi),'ko',...
    'markerfacecolor','k');
xlim([0 2*pi])
ylim([0 2*pi]) 
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
% legend([p1 p2],...
%     'vortices','Wave packet')
% pause(0.05)
%% movie
clf
for i=1:10:length(T)
plot(mod(Y(i,1:M),2*pi),mod(Y(i,M+[1:M]),2*pi),'*','Color',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)],'markerfacecolor',...
    [1-i/numel(T),1-i/numel(T),1-i/numel(T)]);
hold on
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
%% check hamiltonian conservation 
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
P1=zeros(length(T),2);
P2=zeros(length(T),2);
for j=1:length(T)  
    S1=[0 0];    r=[0 0];
    for p=1:length(A)
    S1=S1+A(p).*...
    [Y(j,2*M+p),Y(j,3*M+p)];
    end
    for k=1:length(G)
%         for n=-L:L
%             for m=-L:L
%       r=r+G(k)*[Y(j,4*M+N+k)+2*pi*m,-Y(j,4*M+k)+2*pi*n];
%             end
%         end
      r=r+G(k)*[Y(j,4*M+N+k),-Y(j,4*M+k)];
    end 
    P1(j,:)=S1+r;
%     P2(j,:)=0*S1+r;
end
%
clf
% subplot(2,1,1)
plot(T,P1(:,1),'r') 
% hold on
% subplot(2,1,2)
% hold on
% plot(T,P(:,2),'b')
% plot(T,P2(:,1),'--r')
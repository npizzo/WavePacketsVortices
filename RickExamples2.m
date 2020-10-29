%% correcting the momentum
% N = 15;
cd /Users/npizzo/Documents/Research/Wave_vortex
tstart=cputime;
% N = 1; 
% M = 1;
L = 3;
% put initial points on a lattice
N=10;
M=1; 
init = zeros(4*M+2*N,1);
% we can choose random initial conditions with
% the additional caveat that points aren't too close
[X, Y, ~] = GetPointsRandom(N+M, 2*pi, 2*pi, .25);
for i=1:M
init(i,1) = X(i);
init(M+i,1) = Y(i); 

end
for i=1:N
init(4*M+i,1) = X(i+M);
init(4*M+N+i,1) = Y(i+M);
end
% ----------------
% put initial points on a lattice--------
% No=1+floor(sqrt(No));
% M*No-(M+N);
% ord=randsample(N+M,N+M);
% x=[2*pi/No:2*pi/No:2*pi];
% y=[2*pi/M:2*pi/M:2*pi];
% [X,Y]=meshgrid(x,y);
% rx=reshape(X,[No*M 1]);
% ry=reshape(Y,[No*M 1]);
% %
% init = zeros(4*M+2*N,1);
% for i=1:M
% init(i,1) = mod(rx(ord(i),1),2*pi);
% init(M+i,1) = mod(ry(ord(i),1),2*pi);
% end
% for i=1:N
% init(4*M+i,1) = mod(rx(ord(M+i),1), 2*pi);
% init(4*M+N+i,1) = mod(ry(ord(i+M),1), 2*pi);
% end
% -----------------------
% choose random points-----------
% init(M+1:2*M,1) = 2*pi*rand(M,1);
% init(2*M+1:3*M,1) = 1*ones(length(M),1);
% init(3*M+1:4*M,1) = zeros(M,1);
% init(4*M+1:4*M+N,1) = 2*pi*rand(N,1);
% init(4*M+N+1:4*M+2*N,1) = 2*pi*rand(N,1);
% init(4*M+1:4*M+N,1) = pi+[1:N]/50;
% init(4*M+N+1:4*M+2*N,1) = pi;
% end choose random ---------
init(2*M+1:3*M,1) = 1*ones(length(M),1);
init(3*M+1:4*M,1) = zeros(M,1);
A = 1e-1*ones(M,1); 
G = 1e-1*ones(N,1);
% G=1*( 2*rem(1:N,2) - 1);
g = 1; 
options = odeset('RelTol',1e-7,'AbsTol',1e-7);
[T,Y] = ode15s(@(t,y)WMVN_IP_2( t, y, N, M, G, A, g, L), [0:.1:1000], init, options);
tend=cputime-tstart
%% display data
clf 
% figure
plot(mod(Y(:,1:M), 2*pi), mod(Y(:,M+[1:M]), 2*pi), 'bd',...
    'markerfacecolor','b');
hold on 
plot(mod(Y(:,4*M+[1:N]), 2*pi), mod(Y(:,4*M+N+[1:N]), 2*pi),'k.',...
    'markerfacecolor','k');
xlim([0 2*pi])
ylim([0 2*pi]) 
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
% legend([p1 p2],...
%     'vortices','Wave packet')
% pause(0.05)
%% histogram 
% 
Vx = zeros(length(T)*N,1);
Vy = Vx;
Lt = length(T);
for i = 1 : N
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
%% wavenumber dependence
clf
plot(T,sqrt(Y(:,2*M+1).^2+Y(:,3*M+1).^2))
set(gca,'fontsize',22)
xlabel('t','interpreter','latex')
ylabel('$k$','interpreter','latex')
%% movie
clf
for i=1
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
    ((mod(Y(j,p), 2*pi) - ...
    mod(Y(j,4*M+k), 2*pi)+2*pi*n)*Y(j,3*M+p)-...
    (mod(Y(j,M+p), 2*pi) -...
    mod(Y(j,4*M+N+k), 2*pi)+2*pi*m)*Y(j,2*M+p)))/...
    ((mod(Y(j,p), 2*pi)-mod(Y(j,4*M+k),2*pi)+2*pi*n).^2+...
    (mod(Y(j,p+M),2*pi)-mod(Y(j,4*M+N+k),2*pi) +2*pi*m).^2);
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
      r=r+G(k)*G(l)*log(((mod(Y(j,4*M+k),2*pi)-...
          mod(Y(j,4*M+l),2*pi)+2*pi*n).^2+...
          (mod(Y(j,4*M+N+k), 2*pi) - ...
          mod(Y(j,4*M+N+l),2*pi)+2*pi*m).^2));
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
P = zeros(length(T),2);
for j = 1:length(T)  
    S1=[0 0];    r=[0 0];
    for p=1:length(A)
    S1=S1+A(p).*...
    [Y(j,2*M+p),Y(j,3*M+p)];
    end
    for k=1:length(G)
%         for n=-L:L
%             for m=-L:L
%       r=r+G(k)*[mod(Y(j,4*M+N+k), 2*pi)+2*pi*m,...
%           -mod(Y(j,4*M+k),2*pi)+2*pi*n];
%             end
%         end
      r=r+G(k)*[Y(j,4*M+N+k), - Y(j,4*M+k)];
    end 
    P(j,:)=S1+r;
end
clf

plot(T,P(:,1)/P(1,1),'r') 
hold on
plot(T,P(:,2)/P(1,2),'k')

%% scratch on grid

[X, Y, D] = GetPointsRandom(N+M, 2*pi, 2*pi, 1/2)
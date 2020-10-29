function dpsidx=dPsiIPdx(A,X,Y,x,y,k,l,L)
s=0;
for i=1:length(A)
    for n=-L:L
        for m=-L:L
%     if (x-X(i)+2*pi*n).^2+(y-Y(i)+2*pi*m).^2==0
if (x-X(i)-2*pi*n).^2+(y-Y(i)-2*pi*m).^2==0
        kk = 0;
        ll = 1;
    else 
        kk = 1;
        ll = 0;
end 
%         s=s+kk*A(i)*(2*k(i)*(x-X(i)-2*pi*n).*(y-Y(i)-2*pi*m)-l(i)*...
%         ((x-X(i)-2*pi*n).^2-(y-Y(i)-2*pi*m).^2))./...
%     ((ll+(x-X(i)-2*pi*n).^2+(y-Y(i)-2*pi*m).^2).^2);
    s=s+kk*A(i)*(2*k(i)*(x-X(i)-2*pi*n).*(y-Y(i)-2*pi*m)+l(i)*...
        ((x-X(i)-2*pi*n-y+Y(i)+2*pi*m).*(2*(m+n)*pi-x+X(i)-y+Y(i))))./...
    ((ll+(x-X(i)-2*pi*n).^2+(y-Y(i)-2*pi*m).^2).^2);
        end
    end 
end

dpsidx=1/2/pi*s;
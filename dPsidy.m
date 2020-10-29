function dpsidy=dPsidy(A,X,Y,x,y,k,l)
s=0;
for i=1:length(A)
        if (x-X(i)).^2+(y-Y(i)).^2==0
        kk = 0;
        ll = 1;
    else 
        kk = 1;
        ll = 0;
    end 
    s=s+kk*A(i)*(2*l(i)*(x-X(i)).*(y-Y(i))+k(i)*...
        ((x-X(i)).^2-(y-Y(i)).^2))./...
    ((ll+(x-X(i)).^2+(y-Y(i)).^2).^2);
end

dpsidy=-1/2/pi*s;

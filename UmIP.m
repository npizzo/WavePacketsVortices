function u=UmIP(G,x,y,X,Y,L)
s=0;
for i=1:length(G)
    for n=-L:L
        for m=-L:L
    if (x(i)-X-2*pi*n).^2+(y(i)-Y-2*pi*m).^2==0
        k = 0;
        l = 1;
    else 
        k = 1;
        l = 0;
    end 
    s=s+k*G(i).*((y(i)-Y-2*pi*m)./(l+(x(i)-X-2*pi*n).^2+...
        (y(i)-Y-2*pi*m).^2));
        end
    end
end
u=1/2/pi*s;
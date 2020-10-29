function ud=Ud(G,x,y,X,Y)
s=0;
for i=1:length(G)
    s=s+G(i).*((Y-y(i))./((x(i)-X).^2+(y(i)-Y).^2));
end

ud=1/2/pi*s;
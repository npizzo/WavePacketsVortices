function v=Vm(G,x,y,X,Y)
s=0;
for i=1:length(G)
    if (x(i)-X).^2+(y(i)-Y).^2==0
        k = 0;
        l = 1;
    else 
        k = 1;
        l = 0;
    end 
    s=s+k*G(i).*((X-x(i))./(l+(x(i)-X).^2+(y(i)-Y).^2));
end
v=1/2/pi*s;
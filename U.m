function u=Um(G,x,y,X,Y)
u=1/2/pi*sum(G.*((Y-y)./((x-X).^2+(y-Y).^2)));
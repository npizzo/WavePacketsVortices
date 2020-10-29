function vd=Vd(G,x,y,X,Y)
vd=1/2/pi*sum(G.*((x-X)./((x-X).^2+(y-Y).^2)));
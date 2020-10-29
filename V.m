function v=V(G,x,y,X,Y)
v=1/2/pi*sum(G.*((x-X)./((x-X).^2+(y-Y).^2)));
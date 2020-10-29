%1W1V
function dydt=W1V1(t,y)
dydt=zeros(4,1);
dydt(1)=(-2*y(1)*y(4)*y(3)+y(2)*(-y(4)^2+y(3)^2))/(y(3)^2+y(4)^2)^2;

dydt(2)=(2*y(2)*y(4)*y(3)+y(1)*(-y(4)^2+y(3)^2))/(y(3)^2+y(4)^2)^2;

dydt(3)=y(1)/(2*(y(1)^2+y(2)^2)^(3/4))+...
    2*y(4)*(y(1)*y(4)-y(2)*y(3))/(y(3)^2+y(4)^2)^2-...
    (y(1)+y(4))/(y(3)^2+y(4)^2);

dydt(4)=y(2)/(2*(y(1)^2+y(2)^2)^(3/4))-...
    2*y(4)*(y(2)*y(4)+y(1)*y(3))/(y(3)^2+y(4)^2)^2+...
    (y(2)+y(3))/(y(3)^2+y(4)^2);
a=0:2*pi/100:2*pi;
x1=cos(a);y1=sin(a);
x2=0.01*cos(a);y2=0.01*sin(a);
plot(x1,y1,'r')
axis([-1.5 1.5 -1.5 1.5])
hold on
r=0.25*(1.2+0.25.*cos(3.*a));
x4=-0.0+r.*cos(a);y4=0.0+r.*sin(a);
plot(x4,y4,'r-')

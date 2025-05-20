function og=OriginalBoundaryValue(beta0,M,n,loc)

h=2*pi/n;n1=20;n2=30;
x0=loc(1,1);y0=loc(1,2);
for j=1:n+1
    theta1=(j-1)*h;
    z1=cos(theta1);z2=sin(theta1);
h1=1/n1;h2=2*pi/n2;
[v,theta]=meshgrid(h1/2:h1:1-h1/2,h2/2:h2:2*pi-h2/2);
q=zeros(n2,n1);
for k=1:2*M+1
    if k<=M+1
    q=q+beta0(k)*cos((k-1)*theta);
    else
        q=q+beta0(k)*sin((k-M-1)*theta);
    end
end
x1=x0+v.*q.*cos(theta); x2=y0+v.*q.*sin(theta);
g(j,1)=sum(sum(1./(2.*pi).*(1-x1.^2-x2.^2)./((z1-x1).^2+(z2-x2).^2).*v.*q.^2,2).*h1).*h2;
end
g1=g;

og=g1;
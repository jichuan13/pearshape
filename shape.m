clc
clear
tic,
n=30;M=3; ro=0.3032;
loc=[0.0004,-0.0008];
g=BoundaryValue(n,0,0);
% rand('seed',0)
g=(1+0.01*(2*rand(length(g),1)-1)).*g;
dth=2*pi/n; th=[0:dth:2*pi];
for jj=1:n+1 
   for kk=0:M
       tri(jj,kk+1)=cos(kk*th(jj));
      dtri(jj,kk+1)=-kk*sin(kk*th(jj));
    end
    for kk=1:M
       tri(jj,kk+M+1)=sin(kk*th(jj));
        dtri(jj,kk+M+1)=kk*cos(kk*th(jj));
    end
end
beta=zeros(2*M+1,1); 
beta(1)=ro; 
s=tri*beta;
ds=dtri*beta;


tol=norm(g);kk=1;og1=g;
while(kk<=24820&&tol>=1e-6)
    kk;
beta0=beta;
g1=og1;
og=OriginalBoundaryValue(beta0,M,n,loc);
g2=og;
og1=g2;
GG=g-og;

x0=loc(1,1);y0=loc(1,2);
r=s';
rd=ds';
y1=x0+r.*cos(th);y2=y0+r.*sin(th);
y1d=rd.*cos(th)-r.*sin(th);y2d=rd.*sin(th)+r.*cos(th);
z=cos(th)'+1i*sin(th)';x=y1+1i*y2;  
for k=1:length(th); 
    G=(log(abs(z-x(k)))-log(abs(z./abs(z)-abs(z)*x(k))))/(2*pi);
    GGG1(:,k)=G;
end
z=(1+1.e-7)*z;
for k=1:length(th); 
    G=(log(abs(z-x(k)))-log(abs(z./abs(z)-abs(z)*x(k))))/(2*pi);
    GGG2(:,k)=G; 
end
F=(GGG2-GGG1)/1.e-7;
gg=F*diag(y2d.*cos(th)-y1d.*sin(th))*tri*dth;
dg=GG'*gg*dth;

J=gg;R=-GG;
%阻尼最小二乘法
% t0=svd(inv(J'*J));
% lambda=(1/max(t0)+min(t0))/2;
% beta1=beta0'-(inv(J'*J+lambda*eye(length(beta0)))*(J'*R))';
%高斯牛顿法
% beta1=beta0'-(inv(J'*J)*(J'*R))';
%最速下降法
alpha=0.6;
beta1=beta0'-alpha*(J'*R)';


beta=beta1';

beta
tol=norm(g2-g1)
kk=kk+1
s=tri*beta;
ds=dtri*beta;
end
kk;
a=0:2*pi/n:2*pi;
x0=loc(1,1);y0=loc(1,2);
r=s';
x1=x0+r.*cos(a);y1=y0+r.*sin(a);



hold on 
plot(x1,y1,'b')
domain
toc

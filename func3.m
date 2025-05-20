function  fun=func3(theta,r,c,x0,y0)
x1=x0+r.*cos(theta); x2=y0+r.*sin(theta);
z1=c(1,1);z2=c(1,2);
if (x1==0 & x2==0)
 fun=1./(2.*pi).*r;
else
    fun=1./(2.*pi).*(1-x1.^2-x2.^2)./((z1-x1).^2+(z2-x2).^2).*r;
end
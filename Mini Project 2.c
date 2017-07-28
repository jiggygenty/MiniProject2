#include<stdio.h>
#include<math.h>
#define LSIGMA 10
#define BETA 8/3
#define ROH 99.6
double dx(double y,double x)
{
    return(LSIGMA*(y-x));
}
double dy(double x,double z,double y)
{
    return(x*(ROH-z)-y);
}
double dz(double x,double y,double z)
{
    return(x*y-BETA*z);
}
double nx(double dt,double x,double y,double z)
{
    double xstar,xstar2,xstar3,nx;
    xstar=x+0.5*dt*dx(y,x);
    xstar2=x+0.5*dt*dx(y,xstar);
    xstar3=x+dt*dx(y,xstar2);
    nx=x+(dt/6)*(dx(y,x)+2*dx(y,xstar)+2*dx(y,xstar2)+dx(y,xstar3));
    return(nx);
}
double ny(double dt,double x,double y,double z)
{
    double ystar,ystar2,ystar3,ny;
    ystar=y+0.5*dt*dy(x,z,y);
    ystar2=y+0.5*dt*dy(x,z,ystar);
    ystar3=y+dt*dy(x,z,ystar2);
    ny=y+(dt/6)*(dy(x,z,y)+2*dy(x,z,ystar)+2*dy(x,z,ystar2)+dy(x,z,ystar3));
    return(ny);
}
double nz(double dt,double x,double y,double z)
{
    double zstar,zstar2,zstar3,nz;
    zstar=z+0.5*dt*dz(x,y,z);
    zstar2=z+0.5*dt*dz(x,y,zstar);
    zstar3=z+dt*dz(x,y,zstar2);
    nz=z+(dt/6)*(dz(x,y,z)+2*dz(x,y,zstar)+2*dz(x,y,zstar2)+dz(x,y,zstar3));
    return(nz);
}
void main()
{
FILE *fp;
fp=fopen("lorenzdata.csv","w");
double i,x,y,z,t,tMAX=30,dt=0.01,u,desdel,del,phi;
double sx,sy,sz;
x=1;
y=2;
z=3;

for(t=0;t<=tMAX;t+=dt)
{
    for(i=0;i<2;i++)
    {
    sx=nx(dt,x,y,z);
    sy=ny(dt,x,y,z);
    sz=nz(dt,x,y,z);
    }
fprintf(fp,"%lf,%lf,%lf\n",x,y,z);
printf("%.5lf\t%.5lf\t%.5lf\t%.5lf\n",x,y,z,dt);
x=nx(dt,x,y,z);
y=ny(dt,x,y,z);
z=nz(dt,x,y,z);

if (x>y&&x>z)
    u=x;
else if (y>x&&y>z)
    u=y;
else
    u=z;

if(u==x)
{
desdel=0.0001*(fabs(u)+dt*fabs(dx(y,x)));
del=fabs(sx-u);
phi=desdel/del;
if (phi>1)
dt*=pow(phi,0.2);
if (phi<1)
dt*=pow(phi,0.25);
}
if(u==y)
{
desdel=0.0001*(fabs(u)+dt*fabs(dy(x,z,y)));
del=fabs(sy-u);
phi=desdel/del;
if (phi>1)
dt*=pow(phi,0.2);
if (phi<1)
dt*=pow(phi,0.25);
}
if(u==z)
{
desdel=0.0001*(fabs(u)+dt*fabs(dz(x,y,z)));
del=fabs(sz-u);
phi=desdel/del;
if (phi>1)
dt*=pow(phi,0.2);
if (phi<1)
dt*=pow(phi,0.25);
}
}
fclose(fp);
}

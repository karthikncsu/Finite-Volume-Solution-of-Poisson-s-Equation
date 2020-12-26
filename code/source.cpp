//Calculating the Source term

#include "functions.h"

void source()
{
//-----------------------------------------------------------
//Calculating the source term f(x,y)= exp(-35[(x-x0)^2+(y-y0)^2])
//(x0,y0): coordinates of center of the domain
//-----------------------------------------------------------
double x0=0.0, y0=0.0;
double tot,voltot=0;

for(int j=1;j<jmax; j++)
for(int i=1;i<imax; i++)
{
	x0=x0+xc[i][j]*vol[i][j];
    y0=y0+yc[i][j]*vol[i][j];   
    voltot=voltot+vol[i][j];    
}
x0=x0/voltot;
y0=y0/voltot;

for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
{
	tot=35*((xc[i][j]-x0)*(xc[i][j]-x0)+(yc[i][j]-y0)*(yc[i][j]-y0));
	b[i][j]=exp(-tot)*vol[i][j];
}

}



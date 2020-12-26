//Function for assigning Boundary conditions and source terms

#include "functions.h"

void boundary_conditions()
{
//Right and left boundary condition;
for(int j=1; j<jmax; j++) 
{
if(lb_type==1) usol[0][j]=2*(double)lb_value-usol[1][j];
else usol[0][j]=usol[1][j]-2*(double)lb_value;

if(rb_type==1) usol[imax][j]=2*(double)rb_value-usol[imax-1][j];
else usol[imax][j]=usol[imax-1][j]-2*(double)rb_value;
}

//Top and Bottom boundary condition u=0;
for(int i=1; i<imax; i++) 
{
if(bb_type==1) usol[i][0]=2*bb_value-usol[i][1];
else usol[i][0]=usol[i][1]-2*(double)bb_value;

if(tb_type==1) usol[i][jmax]=2*tb_value-usol[i][jmax-1];
else usol[i][jmax]=usol[i][jmax-1]-2*(double)tb_value;
}

usol[0][0]=0.5*(usol[1][0]+usol[0][1]);
usol[imax][0]=0.5*(usol[imax-1][0]+usol[imax][1]);
usol[imax][jmax]=0.5*(usol[imax-1][jmax]+usol[imax][jmax-1]);
usol[0][jmax]=0.5*(usol[1][jmax]+usol[0][jmax-1]);

}



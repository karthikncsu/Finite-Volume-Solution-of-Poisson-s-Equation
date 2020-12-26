//Program for applying Initial Conditions

#include "functions.h"

void initial_conditions()
{
//----------------------------------------------------------
// Assigning initial conditions
//----------------------------------------------------------

usol=new double* [imax+1];
for (int i=0; i < imax+1; i++) usol[i]= new double[jmax+1];
for(int i=0;i<imax+1;i++) for(int j=0;j<jmax+1;j++) usol[i][j]=0;

b=new double* [imax];
for (int i=0; i < imax; i++) b[i]= new double[jmax];

}


// Function to discretize the poissons equation using Tangent/normal
// Decomposition method

#include "functions.h"

void tangent_normal()
{

source();
boundary_conditions();
	
//----------------------------------------------------------------
//Calculating the coefficients of the K matrix
//----------------------------------------------------------------
a=new double** [imax+1];
for (int i=0; i<imax+1; i++) a[i]= new double* [jmax+1];
for(int i=0; i<imax+1; i++) 
for(int j=0; j<jmax+1; j++)
a[i][j]= new double[5];

for(int i=0;i<imax+1;i++)
for(int j=0;j<jmax+1;j++)
for(int k=0;k<5;k++)
a[i][j][k]=0;

double *bhalf;
bhalf=new double[imax+jmax];

//--------------------------------------
//K Matrix Calculation
//-------------------------------------

//--------------------------------------
//Contribution from average(du) 
//-------------------------------------


for(int i=1;i<imax;i++)
{
double A1[2],A2[2],A3[2],A4[2],Atot[2];
for(int j=1;j<jmax;j++)
{
for(int k=0;k<2;k++)
{	
A1[k]=del[i][j][0][k];
A2[k]=del[i][j][1][k];
A3[k]=del[i-1][j][0][k];
A4[k]=del[i][j-1][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=0.25*(Atot[0]*Atot[0]+Atot[1]*Atot[1])/vol[i][j];
a[i][j][0]=0.25*(-Atot[0]*A3[0]-Atot[1]*A3[1])/vol[i][j];
a[i][j][1]=0.25*(-Atot[0]*A4[0]-Atot[1]*A4[1])/vol[i][j];
a[i][j][3]=0.25*(Atot[0]*A2[0]+Atot[1]*A2[1])/vol[i][j];
a[i][j][4]=0.25*(Atot[0]*A1[0]+Atot[1]*A1[1])/vol[i][j];

for(int k=0;k<2;k++)
{	
A1[k]=del[i+1][j][0][k];
A2[k]=del[i+1][j][1][k];
A3[k]=del[i][j][0][k];
A4[k]=del[i+1][j-1][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(A3[0]*A3[0]+A3[1]*A3[1])/vol[i+1][j];
a[i][j][4]=a[i][j][4]+0.25*(Atot[0]*A3[0]+Atot[1]*A3[1])/vol[i+1][j];

for(int k=0;k<2;k++)
{	
A1[k]=del[i][j+1][0][k];
A2[k]=del[i][j+1][1][k];
A3[k]=del[i-1][j+1][0][k];
A4[k]=del[i][j][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(A4[0]*A4[0]+A4[1]*A4[1])/vol[i][j+1];
a[i][j][3]=a[i][j][3]+0.25*(Atot[0]*A4[0]+Atot[1]*A4[1])/vol[i][j+1];

for(int k=0;k<2;k++)
{	
A1[k]=del[i-1][j][0][k];
A2[k]=del[i-1][j][1][k];
if(i==1) A3[k]=del[0][j][0][k];
else A3[k]=del[i-2][j][0][k];
A4[k]=del[i-1][j-1][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(A1[0]*A1[0]+A1[1]*A1[1])/vol[i-1][j];
a[i][j][0]=a[i][j][0]-0.25*(Atot[0]*A1[0]+Atot[1]*A1[1])/vol[i-1][j];

for(int k=0;k<2;k++)
{	
A1[k]=del[i][j-1][0][k];
A2[k]=del[i][j-1][1][k];
A3[k]=del[i-1][j-1][0][k];
if(j==1) A4[k]=del[i][0][1][k];
else A4[k]=del[i][j-2][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}

a[i][j][2]=a[i][j][2]-0.25*(A2[0]*A2[0]+A2[1]*A2[1])/vol[i][j-1];
a[i][j][1]=a[i][j][1]-0.25*(Atot[0]*A2[0]+Atot[1]*A2[1])/vol[i][j-1];
}
}

//--------------------------------------
//Contribution from (delu.tau)*taudotn*n
//-------------------------------------
for(int i=1;i<imax;i++)
{
double k1[2],k2[2],k3[2],k4[2],ktot[2];	
double A1[2],A2[2],A3[2],A4[2],Atot[2];
double tau[2],dels,taudotA;
for(int j=1;j<jmax;j++)
{
tau[0]=xc[i+1][j]-xc[i][j];
tau[1]=yc[i+1][j]-yc[i][j];
dels=sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
tau[0]=tau[0]/dels;
tau[1]=tau[1]/dels;
taudotA=tau[0]*del[i][j][0][0]+tau[1]*del[i][j][0][1];
for(int k=0;k<2;k++) k1[k]=taudotA*tau[k];

tau[0]=xc[i][j+1]-xc[i][j];
tau[1]=yc[i][j+1]-yc[i][j];
dels=sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
tau[0]=tau[0]/dels;
tau[1]=tau[1]/dels;
taudotA=tau[0]*del[i][j][1][0]+tau[1]*del[i][j][1][1];
for(int k=0;k<2;k++) k2[k]=taudotA*tau[k];

tau[0]=xc[i][j]-xc[i-1][j];
tau[1]=yc[i][j]-yc[i-1][j];
dels=sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
tau[0]=tau[0]/dels;
tau[1]=tau[1]/dels;
taudotA=tau[0]*del[i-1][j][0][0]+tau[1]*del[i-1][j][0][1];
for(int k=0;k<2;k++) k3[k]=taudotA*tau[k];

tau[0]=xc[i][j]-xc[i][j-1];
tau[1]=yc[i][j]-yc[i][j-1];
dels=sqrt(tau[0]*tau[0]+tau[1]*tau[1]);
tau[0]=tau[0]/dels;
tau[1]=tau[1]/dels;
taudotA=+tau[0]*del[i][j-1][1][0]+tau[1]*del[i][j-1][1][1];
for(int k=0;k<2;k++) k4[k]=taudotA*tau[k];

for(int k=0;k<2;k++) ktot[k]=k1[k]+k2[k]-k3[k]-k4[k];
	
for(int k=0;k<2;k++)
{	
A1[k]=del[i][j][0][k];
A2[k]=del[i][j][1][k];
A3[k]=del[i-1][j][0][k];
A4[k]=del[i][j-1][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(Atot[0]*ktot[0]+Atot[1]*ktot[1])/vol[i][j];
a[i][j][0]=a[i][j][0]-0.25*(-ktot[0]*A3[0]-ktot[1]*A3[1])/vol[i][j];
a[i][j][1]=a[i][j][1]-0.25*(-ktot[0]*A4[0]-ktot[1]*A4[1])/vol[i][j];
a[i][j][3]=a[i][j][3]-0.25*(ktot[0]*A2[0]+ktot[1]*A2[1])/vol[i][j];
a[i][j][4]=a[i][j][4]-0.25*(ktot[0]*A1[0]+ktot[1]*A1[1])/vol[i][j];

for(int k=0;k<2;k++)
{	
A1[k]=del[i+1][j][0][k];
A2[k]=del[i+1][j][1][k];
A3[k]=del[i][j][0][k];
A4[k]=del[i+1][j-1][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(-A3[0]*k1[0]-A3[1]*k1[1])/vol[i+1][j];
a[i][j][4]=a[i][j][4]-0.25*(Atot[0]*k1[0]+Atot[1]*k1[1])/vol[i+1][j];

for(int k=0;k<2;k++)
{	
A1[k]=del[i][j+1][0][k];
A2[k]=del[i][j+1][1][k];
A3[k]=del[i-1][j+1][0][k];
A4[k]=del[i][j][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(-A4[0]*k2[0]-A4[1]*k2[1])/vol[i][j+1];
a[i][j][3]=a[i][j][3]-0.25*(Atot[0]*k2[0]+Atot[1]*k2[1])/vol[i][j+1];

for(int k=0;k<2;k++)
{	
A1[k]=del[i-1][j][0][k];
A2[k]=del[i-1][j][1][k];
if(i==1) A3[k]=del[0][j][0][k];
else A3[k]=del[i-2][j][0][k];
A4[k]=del[i-1][j-1][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(-A1[0]*k3[0]-A1[1]*k3[1])/vol[i-1][j];
a[i][j][0]=a[i][j][0]-0.25*(-Atot[0]*k3[0]-Atot[1]*k3[1])/vol[i-1][j];

for(int k=0;k<2;k++)
{	
A1[k]=del[i][j-1][0][k];
A2[k]=del[i][j-1][1][k];
A3[k]=del[i-1][j-1][0][k];
if(j==1) A4[k]=del[i][0][1][k];
else A4[k]=del[i][j-2][1][k];
Atot[k]=A1[k]+A2[k]-A3[k]-A4[k];
}
a[i][j][2]=a[i][j][2]-0.25*(-A2[0]*k4[0]-A2[1]*k4[1])/vol[i][j-1];
a[i][j][1]=a[i][j][1]-0.25*(-Atot[0]*k4[0]-Atot[1]*k4[1])/vol[i][j-1];

}
}

//-----------------------------------
//Contribution from du/dtau*(taudotn)*n
//--------------------------------
for(int j=1;j<jmax; j++)
{
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double taudotn,facearea;
for(int i=0; i<imax; i++)
{
delt1=xc[i+1][j]-xc[i][j];
delt2=yc[i+1][j]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
facearea=(del[i][j][0][0]*del[i][j][0][0]+del[i][j][0][1]+del[i][j][0][1]);
xn1=del[i][j][0][0]/sqrt(facearea);
xn2=del[i][j][0][1]/sqrt(facearea);
taudotn=tau1*xn1+tau2*xn2;
bhalf[i]=-taudotn/dels*sqrt(facearea);
}
for(int i=1;i<imax;i++)
{
a[i][j][2]=a[i][j][2]+bhalf[i]+bhalf[i-1];
a[i][j][0]=a[i][j][0]-bhalf[i-1];
a[i][j][4]=a[i][j][4]-bhalf[i];
}

}

for(int i=1;i<imax; i++)
{
double xc1,xc2,yc1,yc2;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double taudotn,facearea;

for(int j=0; j<jmax; j++)
{
delt1=xc[i][j+1]-xc[i][j];
delt2=yc[i][j+1]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
facearea=(del[i][j][1][0]*del[i][j][1][0]+del[i][j][1][1]+del[i][j][1][1]);
xn1=del[i][j][1][0]/sqrt(facearea);
xn2=del[i][j][1][1]/sqrt(facearea);
taudotn=tau1*xn1+tau2*xn2;
bhalf[j]=-taudotn/dels*sqrt(facearea);
}

for(int j=1;j<jmax;j++)
{
a[i][j][2]=a[i][j][2]+bhalf[j]+bhalf[j-1];
a[i][j][1]=a[i][j][1]-bhalf[j-1];
a[i][j][3]=a[i][j][3]-bhalf[j];

}
}

//------------------------------------------------------------------
//End of K Matrix Calculations
//------------------------------------------------------------------

int Nconv;
double resnorm0=0,resnorm,ratio[Niter];
double **res;
res=new double* [imax];
for(int i=0;i<imax+1;i++) res[i]=new double[jmax];

for(int i=0;i<imax;i++)
for(int j=0;j<jmax;j++)
res[i][j]=0;

double **deltaX;
deltaX=new double* [imax+1];
for(int i=0;i<imax+1;i++) deltaX[i]=new double[jmax+1];

double ***grad;
grad=new double** [imax+1];
for(int i=0;i<imax+1;i++) grad[i]=new double* [jmax+1];
for(int i=0;i<imax+1;i++) for(int j=0;j<jmax+1;j++) grad[i][j]=new double[2];

double D_bar[imax][jmax];

//Algorithm Started
for(int iter=0;iter<Niter;iter++)
{
boundary_conditions();
for(int j=0;j<jmax;j++)
for(int i=0;i<imax;i++)
res[i][j]=0.0;

//---------------------------------------------------------
//Calculating Cell centered gradients using Green's theorem
//---------------------------------------------------------
for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
double f1,f2,f3,f4;
f1 = 0.5*(usol[i][j]+usol[i+1][j]);
f2 = 0.5*(usol[i][j]+usol[i][j+1]);
f3 = 0.5*(usol[i][j]+usol[i-1][j]);
f4 = 0.5*(usol[i][j]+usol[i][j-1]);
for(int k=0;k<2;k++)
{
grad[i][j][k] = f1*del[i][j][0][k] - f3*del[i-1][j][0][k]
			  + f2*del[i][j][1][k] - f4*del[i][j-1][1][k];
grad[i][j][k] = grad[i][j][k]/vol[i][j];
}
}

/*
for(int j=1;j<jmax;j++)
{
double usolgxl,usolgxr;
double f1,f2,f3,f4;
if(lb_type==1) usolgxl=2*lb_value-usol[2][j];
else usolgxl=usol[0][j]-lb_value;
f1 = 0.5*(usol[0][j]+usol[1][j]);
f2 = 0.5*(usol[0][j]+usol[0][j+1]);
f3 = 0.5*(usol[0][j]+usolgxl);
f4 = 0.5*(usol[0][j]+usol[0][j-1]);
for(int k=0;k<2;k++)
{
grad[0][j][k] = f1*del[0][j][0][k] - f3*del[0][j][0][k]
			  + f2*del[0][j][1][k] - f4*del[0][j-1][1][k];
grad[0][j][k] = grad[0][j][k]/vol[0][j];
}
if(rb_type==1) usolgxr=2*rb_value-usol[imax-1][j];
else usolgxr=usol[imax][j]-rb_value;
f1 = 0.5*(usol[imax][j]+usolgxr);
f2 = 0.5*(usol[imax][j]+usol[imax][j+1]);
f3 = 0.5*(usol[imax][j]+usol[imax-1][j]);
f4 = 0.5*(usol[imax][j]+usol[imax][j-1]);
for(int k=0;k<2;k++)
{
grad[imax][j][k] = f1*del[imax][j][0][k] - f3*del[imax-1][j][0][k]
			  + f2*del[imax][j][1][k] - f4*del[imax][j-1][1][k];
grad[imax][j][k] = grad[imax][j][k]/vol[0][j];
}
}

for(int i=1;i<imax;i++)
{
double usolgyb,usolgyt;
double f1,f2,f3,f4;
if(bb_type==1) usolgyb=2*bb_value-usol[i][2];
else usolgyb=usol[i][0]-bb_value;
f1 = 0.5*(usol[i][0]+usol[i+1][0]);
f2 = 0.5*(usol[i][0]+usol[i][1]);
f3 = 0.5*(usol[i][0]+usol[i-1][0]);
f4 = 0.5*(usol[i][0]+usolgyb);
for(int k=0;k<2;k++)
{
grad[i][0][k] = f1*del[i][0][0][k] - f3*del[i-1][0][0][k]
			  + f2*del[i][0][1][k] - f4*del[i-1][0][1][k];
grad[i][0][k] = grad[i][0][k]/vol[i][0];
}
}
*/

//----------------------------------
//Calculating i fluxes
//-----------------------------------
for(int j=1;j<jmax; j++)
{
double *gflux;
gflux=new double[imax];
for(int i=0; i<imax; i++)
{
double f1,f3;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double areaface,taudotn,graddottau;
double gradave1,gradave2;
double gradf1,gradf2;
f1 = usol[i+1][j];
f3 = usol[i][j];
delt1=xc[i+1][j]-xc[i][j];
delt2=yc[i+1][j]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
areaface=sqrt(pow(del[i][j][0][0],2)+pow(del[i][j][0][1],2));
xn1=del[i][j][0][0]/areaface;
xn2=del[i][j][0][1]/areaface;
taudotn=tau1*xn1+tau2*xn2;
gradave1= 0.5*(grad[i][j][0]+grad[i+1][j][0]);
gradave2= 0.5*(grad[i][j][1]+grad[i+1][j][1]);
graddottau = (gradave1*tau1+ gradave2*tau2)*taudotn;
gradf1=(f1-f3)*taudotn*xn1/dels+gradave1-graddottau*xn1;
gradf2=(f1-f3)*taudotn*xn2/dels+gradave2-graddottau*xn2;
gflux[i]=gradf1*del[i][j][0][0]+gradf2*del[i][j][0][1];
}
for(int i=1; i<imax; i++) res[i][j]=res[i][j]+gflux[i]-gflux[i-1];
}

//----------------------------------
//Calculating J fluxes
//-----------------------------------
for(int i=1;i<imax; i++)
{
double *hflux;
hflux=new double[jmax];
for(int j=0;j<jmax; j++)
{
double f2,f4;
double delt1,delt2,dels;
double xn1,xn2;
double tau1,tau2;
double areaface,taudotn,graddottau;
double gradave1,gradave2;
double gradf1,gradf2;
f2 = usol[i][j+1];
f4 = usol[i][j];
delt1=xc[i][j+1]-xc[i][j];
delt2=yc[i][j+1]-yc[i][j];
dels=sqrt(delt1*delt1 + delt2*delt2);
tau1=delt1/dels;                             
tau2=delt2/dels;
areaface=sqrt(pow(del[i][j][1][0],2)+pow(del[i][j][1][1],2));
xn1=del[i][j][1][0]/areaface;
xn2=del[i][j][1][1]/areaface;
taudotn=tau1*xn1+tau2*xn2;
gradave1= 0.5*(grad[i][j][0]+grad[i][j+1][0]);
gradave2= 0.5*(grad[i][j][1]+grad[i][j+1][1]);
graddottau = (gradave1*tau1+ gradave2*tau2)*taudotn;
gradf1=(f2-f4)*taudotn*xn1/dels+gradave1-graddottau*xn1;
gradf2=(f2-f4)*taudotn*xn2/dels+gradave2-graddottau*xn2;
hflux[j]=gradf1*del[i][j][1][0]+gradf2*del[i][j][1][1];
}
for(int j=1; j<jmax; j++) res[i][j]=res[i][j]+hflux[j]-hflux[j-1];
}

resnorm=0;
for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
res[i][j]=-res[i][j]+b[i][j];
resnorm=resnorm+res[i][j]*res[i][j];
}

resnorm=sqrt(resnorm);
if(iter==0) resnorm0=resnorm;	
ratio[iter]=resnorm/resnorm0;
cout<<"Iteration Number: "<<iter<<endl;
cout<<"Absolute Error:"<<resnorm<<endl;
cout<<"Relative Error:"<<ratio[iter]<<endl;
cout<<"-----------------------------------------------------"<<endl;
Nconv=iter+1;
if(ratio[iter]<tol) break;

//------------------------------------------------------------
//Solution using point Jacobian
//------------------------------------------------------------
if(int_method==1)
{
for(int j=1;j<jmax;j++) for(int i=1;i<imax;i++)	deltaX[i][j]=res[i][j]/a[i][j][2];
for(int j=1;j<jmax;j++) for(int i=1;i<imax;i++)	usol[i][j]=deltaX[i][j]+usol[i][j];
}

//------------------------------------------------------------
//Solution using Incomplete LU Decomposition
//------------------------------------------------------------
else if(int_method==2)
{
cout<<"Using Incomplete LU Decomposition"<<endl;
for(int j=1;j<jmax; j++)
for(int i=1;i<imax;i++)
{
	double termi,termj;
	termi=0.0;
	termj=0.0;
	if(i>1) termi=a[i][j][0]*a[i-1][j][4]/D_bar[i-1][j];
	if(j>1) termj=a[i][j][1]*a[i][j-1][3]/D_bar[i][j-1];
	D_bar[i][j]=a[i][j][2]-termi-termj;
}

for(int j=1;j<jmax;j++) for(int i=1;i<imax;i++)	
deltaX[i][j]=(res[i][j]-a[i][j][0]*deltaX[i-1][j]-a[i][j][1]*deltaX[i][j-1])/D_bar[i][j];

for(int j=jmax-1;j>0;j--) for(int i=imax-1;i>0;i--)	
deltaX[i][j]=deltaX[i][j]-(a[i][j][3]*deltaX[i][j+1]+a[i][j][4]*deltaX[i+1][j])/D_bar[i][j];

for(int j=1;j<jmax;j++) for(int i=1;i<imax;i++)	usol[i][j]=deltaX[i][j]+usol[i][j];

for(int j=0;j<jmax+1;j++)
for(int i=0;i<imax+1;i++)deltaX[i][j]=0;
}
//--------------------------------------------
//End of LU Decmposition Algorithm
//--------------------------------------------
}	
//-----------------------------------------
//End of Iteration loop
//-------------------------------------

if(ratio[Nconv-1]<tol)
cout<<"The solver converged"<<endl;
else
{
cout<<"Solution did not converge"<<endl;
cout<<"Maximum number of iterations reached, N="<<Niter<<endl;
}

fstream fout;
string filename;
filename=mesh_name+"_rr";
filename=filename+"_normal_tangent";
if(int_method==1) filename=filename+"_PJ.dat";
else if(int_method==2) filename=filename+"_LU_decomp.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"ZONE I= "<<Nconv<<" F=point"<<endl;
fout<<"Variables = Iteration, Relative_Residue"<<endl;

for(int i=0;i<Nconv;i++)
fout<<i+1<<" "<<log10(ratio[i])<<endl;
fout.close();

filename=mesh_name+"_cputime";
filename=filename+"_normal_tangent";
if(int_method==1) filename=filename+"_PJ.dat";
else if(int_method==2) filename=filename+"_LU_decomp.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"Convergence results for "<<mesh_name;
fout<<" using Normal Tangent Method for spacial discretization";
if(int_method==1) fout<<" and using point jacobi integration method"<<endl;
else if(int_method==2) fout<<" and using incomplete LU decomposition integration method"<<endl;

fout<<"Number of Iterations: "<<Nconv<<endl;
fout.close();
	
}

// Function to discretize the poissons equation using Thin Layer Gradient model

#include "functions.h"

void thin_layer()
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
double D_bar[imax][jmax];

for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
	a[i][j][0]=del[i-1][j][0][0]*del[i-1][j][0][0]+del[i-1][j][0][1]*del[i-1][j][0][1];
	a[i][j][0]=2*a[i][j][0]/(vol[i][j]+vol[i-1][j]);
	
	a[i][j][1]=del[i][j-1][1][0]*del[i][j-1][1][0]+del[i][j-1][1][1]*del[i][j-1][1][1];
	a[i][j][1]=2*a[i][j][1]/(vol[i][j]+vol[i][j-1]);

	a[i][j][3]=del[i][j][1][0]*del[i][j][1][0]+del[i][j][1][1]*del[i][j][1][1];
	a[i][j][3]=2*a[i][j][3]/(vol[i][j]+vol[i][j+1]);
	
	a[i][j][4]=del[i][j][0][0]*del[i][j][0][0]+del[i][j][0][1]*del[i][j][0][1];
	a[i][j][4]=2*a[i][j][4]/(vol[i][j]+vol[i+1][j]);
		
	a[i][j][2]=-a[i][j][0]-a[i][j][1]-a[i][j][3]-a[i][j][4];
}

int Nconv;
double resnorm0=0,resnorm,ratio[Niter];
double **res;
res=new double* [imax];
for(int i=0;i<imax+1;i++) res[i]=new double[jmax];
double **deltaX;
deltaX=new double* [imax+1];
for(int i=0;i<imax+1;i++) deltaX[i]=new double[jmax+1];

//Algorithm Started
for(int iter=0;iter<Niter;iter++)
{
boundary_conditions();
resnorm=0;
for(int j=1; j<jmax; j++)
for(int i=1; i<imax; i++)
{
res[i][j]=-a[i][j][0]*usol[i-1][j]-a[i][j][1]*usol[i][j-1]
		  -a[i][j][3]*usol[i][j+1]-a[i][j][4]*usol[i+1][j]
		  -a[i][j][2]*usol[i][j]+b[i][j];
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
double termi,termj;
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

for(int j=0;j<jmax+1;j++) for(int i=0;i<imax+1;i++) 	deltaX[i][j]=0;

for(int j=1;j<jmax;j++) for(int i=1;i<imax;i++)	
deltaX[i][j]=(res[i][j]-a[i][j][0]*deltaX[i-1][j]-a[i][j][1]*deltaX[i][j-1])/D_bar[i][j];

for(int j=jmax-1;j>0;j--) for(int i=imax-1;i>0;i--)	
deltaX[i][j]=deltaX[i][j]-(a[i][j][3]*deltaX[i][j+1]+a[i][j][4]*deltaX[i+1][j])/D_bar[i][j];

for(int j=1;j<jmax;j++) for(int i=1;i<imax;i++)	usol[i][j]=deltaX[i][j]+usol[i][j];
}
//-------------------------------------------------------
//End of LU Factorization Method
//-------------------------------------------------------
}	

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
filename=filename+"_thin_layer";
if(int_method==1) filename=filename+"_PJ.dat";
else if(int_method==2) filename=filename+"_LU_decomp.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"ZONE I= "<<Nconv<<" F=point"<<endl;
fout<<"Variables = Iteration, Relative_Residue"<<endl;

for(int i=0;i<Nconv;i++)
fout<<i+1<<" "<<log10(ratio[i])<<endl;
fout.close();

filename=mesh_name+"_cputime";
filename=filename+"_thin_layer";
if(int_method==1) filename=filename+"_PJ.dat";
else if(int_method==2) filename=filename+"_LU_decomp.dat";
fout.open(filename.c_str(), ios::out | ios::trunc);
fout<<"Convergence results for "<<mesh_name;
fout<<" using thin layer method for spacial discretization";
if(int_method==1) fout<<" and using point jacobi integration method"<<endl;
else if(int_method==2) fout<<" and using incomplete LU decomposition integration method"<<endl;

fout<<"Number of Iterations: "<<Nconv<<endl;
fout.close();

}

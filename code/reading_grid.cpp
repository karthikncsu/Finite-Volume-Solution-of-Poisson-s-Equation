// Program to read grid and calculating geometry information

#include "functions.h"

void reading_grid(string filename)
{
string line;
ifstream fin;
fin.open(filename.c_str());
    // Check for successful opening of file
    if(fin.fail())
    {
       cout<<"Error Opening File: "<<filename<<endl;
       exit(1);
    }
//----------------------------------------------------------
//Reading the file
//----------------------------------------------------------

getline(fin,line);
fin>>line>>line>>imax>>line>>jmax;
getline(fin,line);

x=new double* [imax+2];
for (int i=0; i < imax+2; i++) x[i]= new double[jmax+2];
y=new double* [imax+2];
for (int i=0; i < imax+2; i++) y[i]= new double[jmax+2];

//------------------------------------------
//Reading the coordinates
//-------------------------------------------
for(int j=1;j<jmax+1;j++) for(int i=1;i<imax+1;i++)	fin>>x[i][j]>>y[i][j];

//Assigning Ghost Cells
for(int j=1; j<jmax+1;j++)
{
x[0][j]=2*x[1][j]-x[2][j];
y[0][j]=2*y[1][j]-y[2][j];
x[imax+1][j]=2*x[imax][j]-x[imax-1][j];
y[imax+1][j]=2*y[imax][j]-y[imax-1][j];
}

for(int i=1; i<imax+1;i++)
{
x[i][0]=2*x[i][1]-x[i][2];
y[i][0]=2*y[i][1]-y[i][2];
x[i][jmax+1]=2*x[i][jmax]-x[i][jmax-1];
y[i][jmax+1]=2*y[i][jmax]-y[i][jmax-1];
}
x[0][0]=2*x[1][0]-x[2][0];
y[0][0]=2*y[0][1]-y[0][2];
x[imax+1][0]=2*x[imax][0]-x[imax-1][0];
y[imax+1][0]=2*y[imax+1][1]-y[imax+1][2];
x[0][jmax+1]=2*x[1][jmax+1]-x[2][jmax+1];
y[0][jmax+1]=2*y[0][jmax]-y[0][jmax-1];
x[imax+1][jmax+1]=2*x[imax][jmax+1]-x[imax-1][jmax+1];
y[imax+1][jmax+1]=2*y[imax+1][jmax]-y[imax+1][jmax-1];

//---------------------------------------------------------
//Calculating the areas of the faces
//---------------------------------------------------------

del=new double*** [imax+1];
for (int i=0; i < imax+1; i++) del[i]= new double** [jmax+1];
for (int i=0; i < imax+1; i++)
for(int j=0;j< jmax+1;j++)
	{
		del[i][j] = new double* [2];
		for(int k=0;k<2;k++) del[i][j][k]= new double[2]; 
	}

for(int j=0;j<jmax+1; j++)
for(int i=0;i<imax+1; i++)
{
del[i][j][0][0]=y[i+1][j+1]-y[i+1][j];
del[i][j][0][1]=-x[i+1][j+1]+x[i+1][j];
del[i][j][1][0]=y[i][j+1]-y[i+1][j+1];
del[i][j][1][1]=-x[i][j+1]+x[i+1][j+1];
}
//----------------------------------------------------------
//Calculating the volumes and center coordinates of the cells
//----------------------------------------------------------
vol=new double* [imax+1];
for (int i=0; i < imax+1; i++) vol[i]= new double[jmax+1];
xc=new double* [imax+1];
for (int i=0; i < imax+1; i++) xc[i]= new double[jmax+1];
yc=new double* [imax+1];
for (int i=0; i < imax+1; i++) yc[i]= new double[jmax+1];
for(int i=0; i<imax+1; i++)
for(int j=0; j<jmax+1; j++)
{
vol[i][j]=triangle_area(x[i+1][j+1],y[i+1][j+1],x[i][j+1],y[i][j+1],x[i][j],y[i][j]);
vol[i][j]=vol[i][j]+triangle_area(x[i+1][j+1],y[i+1][j+1],x[i][j],y[i][j],x[i+1][j],y[i+1][j]);
xc[i][j]=0.25*(x[i][j]+x[i+1][j]+x[i+1][j+1]+x[i][j+1]);
yc[i][j]=0.25*(y[i][j]+y[i+1][j]+y[i+1][j+1]+y[i][j+1]);
}
/*
//-----------------------------------------------------------
//Displaying Geometry information
//-----------------------------------------------------------
cout<<"Coordinates"<<endl;
for(int j=0;j<jmax+1;j++)
for(int i=0;i<imax+1;i++)
cout<<setiosflags(ios::fixed|ios::showpoint)
<<setw(5)<<i<<" "<<setw(5)<<j<<" "
<<setw(10)<<xc[i][j]<<" "<<setw(10)<<yc[i][j]<<endl;
cout<<"-------------------------------------"<<endl;
cout<<"areas"<<endl;
for(int j=1;j<jmax+1;j++)
for(int i=1;i<imax+1;i++)
{
cout<<setiosflags(ios::fixed|ios::showpoint)
<<setw(5)<<i<<" "<<setw(5)<<j<<" "
<<setw(10)<<del[i][j][0][0]<<" "<<setw(10)<<del[i][j][0][1];
cout<<"  "
<<setiosflags(ios::fixed|ios::showpoint)
<<setw(10)<<del[i][j][1][0]<<" "<<setw(10)<<del[i][j][1][1]<<endl;
}

cout<<"-------------------------------------------"<<endl;
cout<<"Volume of each cells"<<endl;
for(int j=0;j<jmax+1;j++)
for(int i=0;i<imax+1;i++)
{
if(vol[i][j]!=0.0004) cout<<i<<" "<<j<<" "<<vol[i][j]-0.0004<<endl;
}
*/


fin.close();
}

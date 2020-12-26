/*Program to create a square grid*/
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>
using namespace std;

int main()
{
int Nx_max,Ny_max;
double lengthx,lengthy,hx,hy,X0,Y0;
cout<<"Enter the points in x direction:";
cin>>Nx_max;
cout<<"Enter the points in y direction:";
cin>>Ny_max;
lengthx=1;
lengthy=1;
X0=0;
Y0=0;
string filename ("squaregrid_");
ostringstream convert;
convert << Nx_max;
filename=filename+convert.str();
filename=filename+'_';
ostringstream convert2;
convert2 << Ny_max;
filename=filename+convert2.str();
filename=filename+".dat";

fstream fout;
fout.open(filename.c_str(), ios::out | ios::trunc);

hx=lengthx/Nx_max;
hy=lengthy/Ny_max;
fout<<"Variables = x, y"<<endl;
fout<<"ZONE I= "<<Nx_max+1<<" J= "<<Ny_max+1<<" F=point"<<endl;
for(int j=0;j<Ny_max+1;j++)
{
	for(int i=0;i<Nx_max+1;i++)
	{
		fout<<X0+i*hx<<" "<<Y0+j*hy<<endl;
	}
}
return(0);
fout.close();
}

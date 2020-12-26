//Function which plots the results

#include "functions.h"

void results()
{
fstream fout,fout1,fout2;
string phifile,gradfile,sourcefile;
phifile=mesh_name+"_phivalues";
gradfile=mesh_name+"_gradient";

if(fv_method==1) phifile=phifile+"_thin_layer";
else if(fv_method==2) phifile=phifile+"_normal_tangent";
if(int_method==1) phifile=phifile+"_PJ.dat";
else if(int_method==2) phifile=phifile+"_LU_decomp.dat";

if(fv_method==1) gradfile=gradfile+"_thin_layer";
else if(fv_method==2) gradfile=gradfile+"_normal_tangent";
if(int_method==1) gradfile=gradfile+"_PJ.dat";
else if(int_method==2) gradfile=gradfile+"_LU_decomp.dat";


fout.open(phifile.c_str(), ios::out | ios::trunc);
fout1.open(gradfile.c_str(), ios::out | ios::trunc);
//fout2.open("sourcefile.dat", ios::out | ios::trunc);

fout<<"Variables = x, y, u"<<endl;
fout<<"ZONE I= "<<imax-1<<" J= "<<jmax-1<<" F=point"<<endl;
fout1<<"Variables = x, y, ugrad, ux, uy"<<endl;
fout1<<"ZONE I= "<<imax-1<<" J= "<<jmax-1<<" F=point"<<endl;
//fout2<<"Variables = x, y, source"<<endl;
//fout2<<"ZONE I= "<<imax-1<<" J= "<<jmax-1<<" F=point"<<endl;

double u1,u2,u3,u4,ux,uy;
for(int j=1;j<jmax;j++)
for(int i=1;i<imax;i++)
{
u1=0.5*(usol[i][j]+usol[i+1][j]);
u2=0.5*(usol[i][j]+usol[i][j+1]);
u3=0.5*(usol[i][j]+usol[i-1][j]);
u4=0.5*(usol[i][j]+usol[i][j-1]);

ux=u1*del[i][j][0][0]+u2*del[i][j][1][0]
  -u3*del[i-1][j][0][0]-u4*del[i][j-1][1][0];
ux=ux/vol[i][j];

uy=u1*del[i][j][0][1]+u2*del[i][j][1][1]
  -u3*del[i-1][j][0][1]-u4*del[i][j-1][1][1];
uy=uy/vol[i][j];

fout<<xc[i][j]<<" "<<yc[i][j]<<" "<<usol[i][j]<<endl;
fout1<<xc[i][j]<<" "<<yc[i][j]<<" "<<sqrt(ux*ux+uy*uy)<<" "<<ux<<" "<<uy<<endl;
//fout2<<xc[i][j]<<" "<<yc[i][j]<<" "<<b[i][j]/vol[i][j]<<endl;

}

fout.close();
fout1.close();
}

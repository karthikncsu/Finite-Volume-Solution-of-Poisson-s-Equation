//Program for solving Poission equation using different Finite Volume Methods

#include "functions.h"
#include <time.h>
#include <sys/time.h>

int imax; // Number of mesh nodes along x direction
int jmax; // Number of mesh nodes along y direction
double ****del; // Node numbers of each element
double **vol; // Coordinates of node points
double **usol; //Solution
double ***a; //Coefficients of K Matrix
double **b; //Source Term
double **x,**y; //X and Y coordinates
double **xc,**yc; //X and Y coordinates of cente of the cell
string mesh_name; //Mesh file name
int fv_method; //Variable indicating finite volume method used
int rb_type; //Boundary condtion type of right boundary
int tb_type; //Boundary condtion type of top boundary
int lb_type; //Boundary condtion type of left boundary
int bb_type; //Boundary condtion type of bottom boundary
double rb_value; //Value of boundary condition for bottom boundary
double tb_value; //Value of boundary condition for bottom boundary
double lb_value; //Value of boundary condition for bottom boundary
double bb_value; //Value of boundary condition for bottom boundary
int int_method; //Integration method
int Niter;
double tol;


int main()
{
double startcputime, endcputime,cpu_time;
startcputime = (double)clock();

read_input();
reading_grid(mesh_name);
initial_conditions();
if(fv_method==1) thin_layer();
else if(fv_method==2) tangent_normal();
cout<<"check2";

 results();

endcputime=(double)clock();
cpu_time=(endcputime-startcputime)/CLOCKS_PER_SEC;

string filename;
fstream fout;
filename=mesh_name+"_cputime";
if(fv_method==1) filename=filename+"_thin_layer";
else if(fv_method==2) filename=filename+"_normal_tangent";
if(int_method==1) filename=filename+"_PJ.dat";
else if(int_method==2) filename=filename+"_LU_decomp.dat";
fout.open(filename.c_str(), ios::out | ios::app);
fout<<"CPU Time: "<<cpu_time<<" sec"<<endl;
fout.close();

return(0);

}

#ifndef HEADER_H
#define HEADER_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <cstdlib>

using namespace std;

extern int imax; // Number of mesh nodes along x direction
extern int jmax; // Number of mesh nodes along y direction
extern double ****del; // Node numbers of each element
extern double **vol; // Coordinates of node points
extern double **usol; //Solution
extern double ***a; //Coefficients of K Matrix
extern double **b; //Source Term
extern double **x,**y; //X and Y coordinates
extern double **xc,**yc; //X and Y coordinates of cente of the cell
extern string mesh_name; //Mesh file name
extern int fv_method; //Variable indicating finite volume method used
extern int rb_type; //Boundary condtion type of right boundary
extern int tb_type; //Boundary condtion type of top boundary
extern int lb_type; //Boundary condtion type of left boundary
extern int bb_type; //Boundary condtion type of bottom boundary
extern double rb_value; //Value of condition for bottom boundary
extern double tb_value; //Value of condition for bottom boundary
extern double lb_value; //Value of condition for bottom boundary
extern double bb_value; //Value of condition for bottom boundary
extern int int_method; //Integration method
extern int Niter;
extern double tol;


void read_input();
void reading_grid(string);
double triangle_area(double,double,double,double,double,double);
void initial_conditions();
void boundary_conditions();
void source();
void thin_layer();
void tangent_normal();
void results();

#endif








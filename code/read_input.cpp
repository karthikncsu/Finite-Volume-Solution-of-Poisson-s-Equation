// Function to read the input file

#include "functions.h"

void read_input()
{
double rb_value_d,lb_value_d,bb_value_d,tb_value_d;
string line;
ifstream fin;
fin.open("input.dat",std::ios::in);

if(fin.fail())
{
cout<<"Error Opening File: "<<"input.dat"<<endl;
exit(1);
}
getline(fin,line);
fin>>mesh_name;
getline(fin,line);
fin>>fv_method;
getline(fin,line);
fin>>rb_type;
getline(fin,line);
fin>>rb_value_d;
rb_value=rb_value_d;
getline(fin,line);
fin>>tb_type;
getline(fin,line);
fin>>tb_value_d;
tb_value=tb_value_d;
getline(fin,line);
fin>>lb_type;
getline(fin,line);
fin>>lb_value_d;
lb_value=lb_value_d;
getline(fin,line);
fin>>bb_type;
getline(fin,line);
fin>>bb_value_d;
bb_value=bb_value_d;
getline(fin,line);
fin>>int_method;
getline(fin,line);
fin>>Niter;
getline(fin,line);
float tol_d;
fin>>tol_d;
tol=tol_d;
getline(fin,line);

cout<<"tolerance: "<<tol<<endl;

fin.close();
}

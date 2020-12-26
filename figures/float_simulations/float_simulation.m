close all
clear all
fileID = fopen('squaregrid_100_100.dat_rr_thin_layer_LU_decomp.dat','r');
formatSpec = '%d %f';
sizeA = [2 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';

x=A(:,1);
y=A(:,2);
fig=figure('Name',['Relative residue for Simulation using float precision']);
title(['Relative residue using float precision']);
xlabel('Number of Iterations')
ylabel('log10(Relative Residue)')
hold on
plot(x,y,'color','k','LineWidth',2);
hold on

clear A,x,y;

fileID1 = fopen('curvilinear.dat_rr_thin_layer_PJ.dat','r');
formatSpec = '%d %f';
sizeA = [2 Inf];
A = fscanf(fileID1,formatSpec,sizeA);
fclose(fileID1);
A=A';
x=A(:,1);
y=A(:,2);
plot(x,y,'--','color','k','LineWidth',2);
hold on

hold on
grid on
hold on

hleg1 = legend('Square Grid 100x100 - LU Decomposition','Curvilinear Grid - Point Jacobi');
set(hleg1,'Location','NorthEast')


close all
clear all
fileID = fopen('curvilinear.dat_rr_thin_layer_PJ.dat','r');
formatSpec = '%d %f';
sizeA = [2 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';

x=A(:,1);
y=A(:,2);
fig=figure('Name',['Relative residue for different methods on Curvilinear Grid']);
title(['Relative Residue for curvilinear grid']);
xlabel('Number of Iterations')
ylabel('log10(Relative Residue)')
hold on
plot(x,y,'color','k','LineWidth',2);
hold on

clear A,x,y;

fileID1 = fopen('curvilinear.dat_rr_thin_layer_LU_decomp.dat','r');
formatSpec = '%d %f';
sizeA = [2 Inf];
A = fscanf(fileID1,formatSpec,sizeA);
fclose(fileID1);
A=A';
x=A(:,1);
y=A(:,2);
plot(x,y,'--','color','k','LineWidth',2);
hold on

clear A,x,y;
fileID3 = fopen('curvilinear.dat_rr_normal_tangent_PJ.dat','r');
formatSpec = '%d %f';
sizeA = [2 Inf];
A = fscanf(fileID3,formatSpec,sizeA);
fclose(fileID3);
A=A';

flag=1;
indx=size(A)/20;
indx1=round(indx(1));

for i=1:size(A)
    if(mod(i,indx1)==0)
    x1(flag)=A(i,1);
    y1(flag)=A(i,2);
    flag=flag+1;
    end   
end
plot(x1,y1,'--o','color','k','LineWidth',2);
hold on

clear A,x,y;

fileID4 = fopen('curvilinear.dat_rr_normal_tangent_LU_decomp.dat','r');
formatSpec = '%d %f';
sizeA = [2 Inf];
A = fscanf(fileID4,formatSpec,sizeA);
fclose(fileID4);
A=A';

flag=1;
indx=size(A)/20;
indx1=round(indx(1));

for i=1:size(A)
    if(mod(i,indx1)==0)
    x1(flag)=A(i,1);
    y1(flag)=A(i,2);
    flag=flag+1;
    end   
end
plot(x1,y1,'*','color','k','LineWidth',2);
hold on
grid on
hold on

hleg1 = legend('Thin Layer-Point Jacobi','Thin Layer - Incomplete LU','Normal/Tangent - Point Jacobi','Normal/Tangent - Incomplete LU');
set(hleg1,'Location','NorthEast')


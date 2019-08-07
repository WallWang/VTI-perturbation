#include<math.h>
#include<stdio.h>
#define nx 901
#define nz 901
#define dx 10
#define dz 10

float nb[nx][nz]={0.0};
float np[nx][nz]={0.0};
float na[nx][nz]={0.0};

int main(){
    int i,j;
    double x,z,R,r;
    FILE *f1,*f2,*f3;
    f1=fopen("./delta_background.dat","w");
    f2=fopen("./delta_perturb.dat","w");
    f3=fopen("./delta_all.dat","w");

//    for(i= 0;i<nx;i++)
//        for(j= 0;j<nz;j++){
//            nb[i][j]= 0.1;//+0.00001*i+0.00001*j;
//        }

    R= 400;

    for(i= 0;i<nx;i++)
        for(j= 0;j<nz;j++)
        {
            x=i*dx;
            z=j*dz;
            r= sqrt(pow(x-4500,2)+pow(z-3000,2));
			nb[i][j]=0.1+0.1*exp(-0.2*pow(r/R,2));
            np[i][j]= 0.*exp(-0.2*pow(r/R,2));
            na[i][j]= nb[i][j]+ np[i][j];
            fwrite(&nb[i][j],sizeof(float),1,f1);
            fwrite(&np[i][j],sizeof(float),1,f2);
            fwrite(&na[i][j],sizeof(float),1,f3);
        }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    return(0);
}

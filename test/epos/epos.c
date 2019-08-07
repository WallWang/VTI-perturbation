#include<math.h>
#include<stdio.h>
#define nx 901
#define nz 901
#define dx 10
#define dz 10

float eb[nx][nz]={0.0};
float ep[nx][nz]={0.0};
float ea[nx][nz]={0.0};

int main(){
    int i,j;
    double x,z,R,r;
    FILE *f1,*f2,*f3;
    f1=fopen("./epos_background.dat","w");
    f2=fopen("./epos_perturb.dat","w");
    f3=fopen("./epos_all.dat","w");

//    for(i= 0;i<nx;i++)
//        for(j= 0;j<nz;j++){
//            eb[i][j]= 0.2;
//        }

    R= 400;

    for(i= 0;i<nx;i++)
        for(j= 0;j<nz;j++)
        {
            x=i*dx;
            z=j*dz;
            r= sqrt(pow(x-4500,2)+pow(z-3000,2));
			eb[i][j]= 0.2+0.2*exp(-0.2*pow(r/R,2));
            ep[i][j]= 0.*exp(-0.2*pow(r/R,2));
            ea[i][j]= eb[i][j]+ ep[i][j];
            fwrite(&eb[i][j],sizeof(float),1,f1);
            fwrite(&ep[i][j],sizeof(float),1,f2);
            fwrite(&ea[i][j],sizeof(float),1,f3);
        }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    return(0);
}

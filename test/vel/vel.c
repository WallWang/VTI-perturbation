#include<math.h>
#include<stdio.h>
#define nx 901
#define nz 901
#define dx 10
#define dz 10

float vb[nx][nz]={0.0};
float vp[nx][nz]={0.0};
float va[nx][nz]={0.0};

int main(){
    int i,j;
    double x,z,R,r;
    FILE *f1,*f2,*f3;
    f1=fopen("./vel_background.dat","w");
    f2=fopen("./vel_perturb.dat","w");
    f3=fopen("./vel_all.dat","w");

//    for(i= 0;i<nx;i++)
//        for(j= 0;j<nz;j++){
//            vb[i][j]= 2000+i+0.3*j;
 //       }

    R= 400;

    for(i= 0;i<nx;i++)
        for(j= 0;j<nz;j++)
        {
            x=i*dx;
            z=j*dz;
            r= sqrt(pow(x-3500,2)+pow(z-3000,2));
			vb[i][j]=3000+0.*exp(-0.2*pow(r/R,2));
            vp[i][j]= 300.*exp(-0.2*pow(r/R,2));
	//	    vp[i][j]=0.0+0.01*x+0.01*z;
            va[i][j]= vb[i][j]+ vp[i][j];
			fwrite(&vb[i][j],sizeof(float),1,f1);
            fwrite(&vp[i][j],sizeof(float),1,f2);
            fwrite(&va[i][j],sizeof(float),1,f3);
        }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    return(0);
}

#include<math.h>
#include<stdio.h>
#define nx 901
#define nz 901

int main(){
    int i,j;
	float bepos,pepos,aepos;
	float benda,penda,aenda;
	float bdelta,pdelta,adelta;
	FILE *fe1,*fe2,*fe3;
	FILE *fn1,*fn2,*fn3;
    FILE *f1,*f2,*f3;
	fe1=fopen("../../../VTI_martix_all_enda/test_stereo2d/epos/epos_background.dat","r");
	fe2=fopen("../../../VTI_martix_all_enda/test_stereo2d/epos/epos_perturb.dat","r");
	fe3=fopen("../../../VTI_martix_all_enda/test_stereo2d/epos/epos_all.dat","r");
	fn1=fopen("../../../VTI_martix_all_enda/test_stereo2d/enda/enda_background.dat","r");
	fn2=fopen("../../../VTI_martix_all_enda/test_stereo2d/enda/enda_perturb.dat","r");
	fn3=fopen("../../../VTI_martix_all_enda/test_stereo2d/enda/enda_all.dat","r");
    f1=fopen("./delta_background.dat","w");
    f2=fopen("./delta_perturb.dat","w");
    f3=fopen("./delta_all.dat","w");
    for(i= 0;i<nx;i++)
        for(j= 0;j<nz;j++)
        {
			fread(&bepos,sizeof(float),1,fe1);
			fread(&pepos,sizeof(float),1,fe2);
			fread(&aepos,sizeof(float),1,fe3);
			fread(&benda,sizeof(float),1,fn1);
			fread(&penda,sizeof(float),1,fn2);
			fread(&aenda,sizeof(float),1,fn3);
			bdelta=bepos-benda;
			pdelta=pepos-penda;
			adelta=aepos-aenda;
			fwrite(&bdelta,sizeof(float),1,f1);
			fwrite(&pdelta,sizeof(float),1,f2);
			fwrite(&adelta,sizeof(float),1,f3);
        }
    fclose(f1);
    fclose(f2);
    fclose(f3);
	fclose(fe1);
	fclose(fe2);
	fclose(fe3);
	fclose(fn1);
	fclose(fn2);
	fclose(fn3);
    return(0);
}

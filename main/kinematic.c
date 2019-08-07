#include "su.h"
#include "stereo_struct.h"
#include "wx.h"
#include "iofile.h"
#include "RaySub.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" 								",
" 								",
" Required parameters:                                          ",
" velocity=                   velocityfield filename        ",
" anisotropic_para_epos=      eposfield filename            ",
" anisotropic_para_delta=     deltafield filename            ",
" raypath_filename=           record raypath to file         ",
" 								",
" Optional Parameters:                                                  ",
" 								",
" fx=0               first lateral sample in velocity          ",
" fz=0               first depth sample in velocity          ",
" dx=                the grid space for velocity model          ",
" dz=                the grid space for velocity model          ",
" nx=                number of lateral samples in velocity          ",
" nz=                number of depth samples in velocity          ",
" ds=                integration step                            ",
" theta=			 initial scattering angle			       	",
" x_start=		 	 lateral location of source  	",
" z_start=			 depth of source			       	",
" 								",
" Notes:                                                                ",
"	stereotomography     					",
NULL};
/**************** end self doc ***********************************/
int main(int argc, char **argv)
{	
	initargs(argc,argv);
	requestdoc(0);

	int i,j;
	int nx,nz;
	float dx,dz;
	float fx,fz;
	float x_start,z_start;
	float theta,ds;
	double vp0,si2,sic;
	double px,pz;
	double vel0,dvdx,dvdz,uxx,uxz,uzz;
	double epos0,dedx,dedz,epxx,epxz,epzz;
	double delta0,deldx,deldz,delxx,delxz,delzz;
	char *raypath_filename;
	char *vel_filename;
	char *epos_filename;
	char *delta_filename;
	BGfield **vel,**epos,**delta;
	Geo2d geo2dv;
	Ray ray;
	Source sc;
	
	//read parameter
	MUSTGETPARSTRING("velocity",&vel_filename);
	MUSTGETPARSTRING("anisotropic_para_epos",&epos_filename);
	MUSTGETPARSTRING("anisotropic_para_delta",&delta_filename);
	MUSTGETPARSTRING("raypath_filename",&raypath_filename);
	if(!getparint("nx",&nx))         nx=1;
	if(!getparint("nz",&nz))         nz=1;
	if(!getparfloat("dx",&dx))       dx=1.0;
	if(!getparfloat("dz",&dz))       dz=1.0;
	if(!getparfloat("fx",&fx))       fx=0;
	if(!getparfloat("fz",&fz))       fz=0;
	if(!getparfloat("ds",&ds))       ds=0.001;
	if(!getparfloat("theta",&theta)) theta=150.;
	if(!getparfloat("x_start",&x_start))  x_start=0.;
	if(!getparfloat("z_start",&z_start))  z_start=0.;
	checkpars();

	//define observation system
	geo2dv.dt=ds;
	geo2dv.nx=nx;  geo2dv.nz=nz;
	geo2dv.fx=fx;  geo2dv.fz=fz;
	geo2dv.dx=dx;  geo2dv.dz=dz;
	geo2dv.xmin=fx;  geo2dv.zmin=fz;
	geo2dv.xmax=fx+(nx-1)*dx;
	geo2dv.zmax=fz+(nz-1)*dz;

	//read vel epos and delta
	FILE *file_vel,*file_epos,*file_delta,*file_path;
	file_vel  = fopen(vel_filename,"r");
	file_epos = fopen(epos_filename,"r");
	file_delta = fopen(delta_filename,"r");
	vel=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	epos=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	delta=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++){
		  fread(&vel[i][j].u,sizeof(float),1,file_vel);
		  fread(&epos[i][j].u,sizeof(float),1,file_epos);
		  fread(&delta[i][j].u,sizeof(float),1,file_delta);
	  }
	fclose(file_vel);
	fclose(file_epos);
	fclose(file_delta);

	//calculate second derivative
	dv2(geo2dv,vel);
	dv2(geo2dv,epos);
	dv2(geo2dv,delta);

	ray.rs=(RayStep*)ealloc1(999999,sizeof(RayStep));
	vel2Interp(geo2dv,vel,x_start,z_start,&vel0,&dvdx,&dvdz,&uxx,&uxz,&uzz);
	vel2Interp(geo2dv,epos,x_start,z_start,&epos0,&dedx,&dedz,&epxx,&epxz,&epzz);
	vel2Interp(geo2dv,delta,x_start,z_start,&delta0,&deldx,&deldz,&delxx,&delxz,&delzz);
	si2=sin(angle2radian(180-theta))*sin(angle2radian(180-theta));
	sic=sin(2*angle2radian(180-theta))*sin(2*angle2radian(180-theta));
	vp0=vel0*sqrt(0.5+epos0*si2+0.5*sqrt(pow(1+2*epos0*si2,2)-2*(epos0-delta0)*sic));//calculate phase velocity
	px=sin(angle2radian(180-theta))/vp0;
	pz=cos(angle2radian(180-theta))/vp0;
	printf("Vphase=%lf,Px0=%e,Pz0=%e\n",vp0,px,pz);
	sc.x=x_start;
	sc.z=z_start;
	sc.px=px;
	sc.pz=pz;
	kinematic_ray_tracing(&ray,sc,geo2dv,vel,epos,delta,99);
	
	file_path=fopen(raypath_filename,"w");
//	for(i=0;i<ray.nrs;i++){
//		printf("------------I=%d---------------\n",i);
//		printf("X=%lf,Z=%lf,Px=%e,Pz=%e,T=%lf\n",ray.rs[i].x,ray.rs[i].z,ray.rs[i].px,ray.rs[i].pz,ray.rs[i].t);
//	}
    printf("==================================================\n");
	printf("T Lp: X=%lf,Z=%lf,Px=%1.6e,Pz=%1.6e,T=%lf\n",ray.rs[ray.nrs-1].x,ray.rs[ray.nrs-1].z,ray.rs[ray.nrs-1].px,ray.rs[ray.nrs-1].pz,ray.rs[ray.nrs-1].t);
	printf("==================================================\n");
	draw_raypath(file_path,ray);
	free(vel);
	free(epos);
	free(delta);
	free(ray.rs);
	fclose(file_path);
	return 0;
}

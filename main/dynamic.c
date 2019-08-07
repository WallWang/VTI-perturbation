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
" velocity=                  velocityfield filename        ",
" anisotropic_para_epos=     eposfield filename            ",
" anisotropic_para_delta=    deltafield filename            ",
" velocity_perturb=          perturbation of velocityfield filename    ",
" anisotropic_para_depos=    perturbation of eposfield filename        ",
" anisotropic_para_ddelta=   perturbation of endafield filename        ",
" center_raypath_filename=   record center raypath to file        ",
" perturb_raypath_filename=  record perturbation raypath to file        ",
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
" delta_theta =		 perturbation of initial scattering angle			       	",
" delta_x =		 	 perturbation of lateral location of source  	",
" delta_z =			 perturbation of depth of source			       	",
" 								",
" Notes:                                                                ",
"	stereotomography     					",
NULL};
/**************** end self doc ***********************************/
int main(int argc, char **argv)
{	
	initargs(argc,argv);
	requestdoc(0);

	int i,j,k;
	int nx,nz;
	float dx,dz;
	float fx,fz;
	float x_start,z_start;
	float delta_x,delta_z;
	float theta,ds;
	float del_theta;
	double Q[4][4];
	double vp0,si2,sic;
	double px0,pz0;
	double vel0,dvdx,dvdz,uxx,uxz,uzz;
	double epos0,dedx,dedz,epxx,epxz,epzz;
	double delta0,deldx,deldz,delxx,delxz,delzz;
	double vpertu,dvpdx,dvpdz,puxx,puxz,puzz;
	double epertu,depdx,depdz,pexx,pexz,pezz;
	double delpertu,delpdx,delpdz,delpxx,delpxz,delpzz;
	double delta_w1,delta_w2,delta_w3,delta_w4;
	double s_xv=0.0,s_zv=0.0,s_pxv=0.0,s_pzv=0.0;
	double t_x,t_z,t_px,t_pz;
	double dhdx,dhdz,dhdpx,dhdpz,dhdt,dtu;
	double bepos,px,pz,px2,pz2,v2,v3,v4;
	char *center_raypath;
	char *perturb_raypath;
	char *vel_filename;
	char *epos_filename;
	char *delta_filename;
	char *vel_pertrub_filename;
	char *epos_perturb_filename;
	char *delta_perturb_filename;
	BGfield **vel,**epos,**delta;
	BGfield **dvel,**depos,**ddel;
	Geo2d geo2dv;
	Ray cray,pray;
	Source sc;
	
	//read parameter
	MUSTGETPARSTRING("velocity",&vel_filename);
	MUSTGETPARSTRING("anisotropic_para_epos",&epos_filename);
	MUSTGETPARSTRING("anisotropic_para_delta",&delta_filename);
	MUSTGETPARSTRING("velocity_perturb",&vel_pertrub_filename);
	MUSTGETPARSTRING("anisotropic_para_depos",&epos_perturb_filename);
	MUSTGETPARSTRING("anisotropic_para_ddelta",&delta_perturb_filename);
	MUSTGETPARSTRING("center_raypath_filename",&center_raypath);
	MUSTGETPARSTRING("perturb_raypath_filename",&perturb_raypath);
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
	if(!getparfloat("delta_theta",&del_theta)) del_theta=0.;
	if(!getparfloat("delta_x",&delta_x))  delta_x=0.;
	if(!getparfloat("delta_z",&delta_z))  delta_z=0.;
	checkpars();

	geo2dv.dt=ds;
	geo2dv.nx=nx;  geo2dv.nz=nz;
	geo2dv.fx=fx;  geo2dv.fz=fz;
	geo2dv.dx=dx;  geo2dv.dz=dz;
	geo2dv.xmin=fx;  geo2dv.zmin=fz;
	geo2dv.xmax=fx+(nx-1)*dx;
	geo2dv.zmax=fz+(nz-1)*dz;
	del_theta=angle2radian(del_theta);
	//read vel epos and delta
	FILE *file_vel,*file_epos,*file_delta,*file_dv,*file_de,*file_dd;
	file_vel  = fopen(vel_filename,"r");
	file_epos = fopen(epos_filename,"r");
	file_delta = fopen(delta_filename,"r");
	file_dv   = fopen(vel_pertrub_filename,"r");
	file_de   = fopen(epos_perturb_filename,"r");
	file_dd   = fopen(delta_perturb_filename,"r");
	vel=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	epos=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	delta=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	dvel=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	depos=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	ddel=(BGfield**)ealloc2(nz,nx,sizeof(BGfield));
	for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++){
		  fread(&vel[i][j].u,sizeof(float),1,file_vel);
		  fread(&epos[i][j].u,sizeof(float),1,file_epos);
		  fread(&delta[i][j].u,sizeof(float),1,file_delta);
		  fread(&dvel[i][j].u,sizeof(float),1,file_dv);
		  fread(&depos[i][j].u,sizeof(float),1,file_de);
		  fread(&ddel[i][j].u,sizeof(float),1,file_dd);
	  }
	fclose(file_vel);
	fclose(file_epos);
	fclose(file_delta);
	fclose(file_dv);
	fclose(file_de);
	fclose(file_dd);
	
	dv2(geo2dv,vel);
	dv2(geo2dv,epos);
	dv2(geo2dv,delta);
	dv2(geo2dv,dvel);
	dv2(geo2dv,depos);
	dv2(geo2dv,ddel);
	
	cray.rs=(RayStep*)ealloc1(999999,sizeof(RayStep));
	pray.rs=(RayStep*)ealloc1(999999,sizeof(RayStep));
	vel2Interp(geo2dv,vel,x_start,z_start,&vel0,&dvdx,&dvdz,&uxx,&uxz,&uzz);
	vel2Interp(geo2dv,epos,x_start,z_start,&epos0,&dedx,&dedz,&epxx,&epxz,&epzz);
	vel2Interp(geo2dv,delta,x_start,z_start,&delta0,&deldx,&deldz,&delxx,&delxz,&delzz);
	si2=sin(angle2radian(180-theta))*sin(angle2radian(180-theta));
	sic=sin(2*angle2radian(180-theta))*sin(2*angle2radian(180-theta));
	vp0=vel0*sqrt(0.5+epos0*si2+0.5*sqrt(pow(1+2*epos0*si2,2)-2*(epos0-delta0)*sic));
	px0=sin(angle2radian(180-theta))/vp0;
	pz0=cos(angle2radian(180-theta))/vp0;
	printf("Vphase=%lf,Px0=%e,Pz0=%e\n",vp0,px0,pz0);

	sc.x=x_start;
	sc.z=z_start;
	sc.px=px0;
	sc.pz=pz0;
	sc.t=0.0;
	dynamic_ray_tracing(&cray,sc,geo2dv,vel,epos,delta,99);
	pray.nrs=0;

	for(i=0;i<cray.nrs-1;i++){
		vel2Interp(geo2dv,vel,cray.rs[i].x,cray.rs[i].z,&vel0,&dvdx,&dvdz,&uxx,&uxz,&uzz);
		vel2Interp(geo2dv,epos,cray.rs[i].x,cray.rs[i].z,&epos0,&dedx,&dedz,&epxx,&epxz,&epzz);
		vel2Interp(geo2dv,delta,cray.rs[i].x,cray.rs[i].z,&delta0,&deldx,&deldz,&delxx,&delxz,&delzz);
		vel2Interp(geo2dv,dvel,cray.rs[i].x,cray.rs[i].z,&vpertu,&dvpdx,&dvpdz,&puxx,&puxz,&puzz);
		vel2Interp(geo2dv,depos,cray.rs[i].x,cray.rs[i].z,&epertu,&depdx,&depdz,&pexx,&pexz,&pezz);
		vel2Interp(geo2dv,ddel,cray.rs[i].x,cray.rs[i].z,&delpertu,&delpdx,&delpdz,&delpxx,&delpxz,&delpzz);
		px=cray.rs[i].px;
		pz=cray.rs[i].pz;
		v2=vel0*vel0;
		v3=v2*vel0;
		v4=v2*v2;
		px2=px*px;
		pz2=pz*pz;
		bepos=1+2*epos0;
		
		delta_w1=(2*vel0*px*(bepos-4*(epos0-delta0)*v2*pz2))*vpertu+2*v2*px*(1-v2*pz2)*epertu+2*v4*px*pz2*delpertu;
		delta_w2=(2*vel0*pz*(1-4*(epos0-delta0)*v2*px2))*vpertu-2*v4*px2*pz*epertu+2*v4*px2*pz*delpertu;
		delta_w3=(4*v2*(epos0-delta0)*px2*pz2-bepos*px2-pz2)*vel0*dvpdx-((bepos*px2+pz2-12*v2*(epos0-delta0)*px2*pz2)*dvdx+2*vel0*px2*dedx-4*v3*px2*pz2*(dedx-deldx))*vpertu+2*vel0*px2*(2*v2*pz2-1)*dvdx*epertu+v2*px2*(v2*pz2-1)*depdx-4*v3*px2*pz2*dvdx*delpertu-v4*px2*pz2*delpdx;
		delta_w4=(4*v2*(epos0-delta0)*px2*pz2-bepos*px2-pz2)*vel0*dvpdz-((bepos*px2+pz2-12*v2*(epos0-delta0)*px2*pz2)*dvdz+2*vel0*px2*dedz-4*v3*px2*pz2*(dedz-deldz))*vpertu+2*vel0*px2*(2*v2*pz2-1)*dvdz*epertu+v2*px2*(v2*pz2-1)*depdz-4*v3*px2*pz2*dvdz*delpertu-v4*px2*pz2*delpdz;

		t_x =cray.rs[i].Q[0][0]*(delta_x+s_xv)+cray.rs[i].Q[0][1]*(delta_z+s_zv)+cray.rs[i].Q[0][2]*(-pz0*del_theta+s_pxv)+cray.rs[i].Q[0][3]*(px0*del_theta+s_pzv);
		t_z =cray.rs[i].Q[1][0]*(delta_x+s_xv)+cray.rs[i].Q[1][1]*(delta_z+s_zv)+cray.rs[i].Q[1][2]*(-pz0*del_theta+s_pxv)+cray.rs[i].Q[1][3]*(px0*del_theta+s_pzv);
		t_px=cray.rs[i].Q[2][0]*(delta_x+s_xv)+cray.rs[i].Q[2][1]*(delta_z+s_zv)+cray.rs[i].Q[2][2]*(-pz0*del_theta+s_pxv)+cray.rs[i].Q[2][3]*(px0*del_theta+s_pzv);
		t_pz=cray.rs[i].Q[3][0]*(delta_x+s_xv)+cray.rs[i].Q[3][1]*(delta_z+s_zv)+cray.rs[i].Q[3][2]*(-pz0*del_theta+s_pxv)+cray.rs[i].Q[3][3]*(px0*del_theta+s_pzv);
	
		s_xv  += (+cray.rs[i+1].Q[2][2]*delta_w1+cray.rs[i+1].Q[3][2]*delta_w2-cray.rs[i+1].Q[0][2]*delta_w3-cray.rs[i+1].Q[1][2]*delta_w4)*(cray.rs[i+1].s-cray.rs[i].s);
		s_zv  += (+cray.rs[i+1].Q[2][3]*delta_w1+cray.rs[i+1].Q[3][3]*delta_w2-cray.rs[i+1].Q[0][3]*delta_w3-cray.rs[i+1].Q[1][3]*delta_w4)*(cray.rs[i+1].s-cray.rs[i].s);
		s_pxv += (-cray.rs[i+1].Q[2][0]*delta_w1-cray.rs[i+1].Q[3][0]*delta_w2+cray.rs[i+1].Q[0][0]*delta_w3+cray.rs[i+1].Q[1][0]*delta_w4)*(cray.rs[i+1].s-cray.rs[i].s);
		s_pzv += (-cray.rs[i+1].Q[2][1]*delta_w1-cray.rs[i+1].Q[3][1]*delta_w2+cray.rs[i+1].Q[0][1]*delta_w3+cray.rs[i+1].Q[1][1]*delta_w4)*(cray.rs[i+1].s-cray.rs[i].s);

		if(cray.rs[i].x+t_x<=geo2dv.xmax&&cray.rs[i].x+t_x>=geo2dv.xmin&&cray.rs[i].z+t_z<=geo2dv.zmax&&cray.rs[i].z+t_z>=geo2dv.zmin){
			pray.rs[pray.nrs].x = cray.rs[i].x+t_x;
			pray.rs[pray.nrs].z = cray.rs[i].z+t_z;
			pray.rs[pray.nrs].px= cray.rs[i].px+t_px;
			pray.rs[pray.nrs].pz= cray.rs[i].pz+t_pz;
			pray.rs[pray.nrs].t = cray.rs[i].t;
			pray.nrs++;
		}
	}
	vel2Interp(geo2dv,vel,cray.rs[cray.nrs-1].x,cray.rs[cray.nrs-1].z,&vel0,&dvdx,&dvdz,&uxx,&uxz,&uzz);
	vel2Interp(geo2dv,epos,cray.rs[cray.nrs-1].x,cray.rs[cray.nrs-1].z,&epos0,&dedx,&dedz,&epxx,&epxz,&epzz);
	vel2Interp(geo2dv,delta,cray.rs[cray.nrs-1].x,cray.rs[cray.nrs-1].z,&delta0,&deldx,&deldz,&delxx,&delxz,&delzz);
	px=cray.rs[cray.nrs-1].px;
	pz=cray.rs[cray.nrs-1].pz;
	v2=vel0*vel0;
	v3=v2*vel0;
	v4=v2*v2;
	px2=px*px;
	pz2=pz*pz;
	bepos=1+2*epos0;

	//extra raytracing
	t_x =cray.rs[cray.nrs-1].Q[0][0]*(delta_x+s_xv)+cray.rs[cray.nrs-1].Q[0][1]*(delta_z+s_zv)+cray.rs[cray.nrs-1].Q[0][2]*(-pz0*del_theta+s_pxv)+cray.rs[cray.nrs-1].Q[0][3]*(px0*del_theta+s_pzv);
	t_z =cray.rs[cray.nrs-1].Q[1][0]*(delta_x+s_xv)+cray.rs[cray.nrs-1].Q[1][1]*(delta_z+s_zv)+cray.rs[cray.nrs-1].Q[1][2]*(-pz0*del_theta+s_pxv)+cray.rs[cray.nrs-1].Q[1][3]*(px0*del_theta+s_pzv);
	t_px=cray.rs[cray.nrs-1].Q[2][0]*(delta_x+s_xv)+cray.rs[cray.nrs-1].Q[2][1]*(delta_z+s_zv)+cray.rs[cray.nrs-1].Q[2][2]*(-pz0*del_theta+s_pxv)+cray.rs[cray.nrs-1].Q[2][3]*(px0*del_theta+s_pzv);
	t_pz=cray.rs[cray.nrs-1].Q[3][0]*(delta_x+s_xv)+cray.rs[cray.nrs-1].Q[3][1]*(delta_z+s_zv)+cray.rs[cray.nrs-1].Q[3][2]*(-pz0*del_theta+s_pxv)+cray.rs[cray.nrs-1].Q[3][3]*(px0*del_theta+s_pzv);

	dhdx=v2*px*(bepos-2*v2*pz2*(epos0-delta0));
	dhdz=v2*pz*(1-2*v2*px2*(epos0-delta0));
	dhdpx=(4*v2*px2*pz2*(epos0-delta0)-bepos*px2-pz2)*vel0*dvdx-v2*px2*dedx+v4*px2*pz2*(dedx-deldx);
	dhdpz=(4*v2*px2*pz2*(epos0-delta0)-bepos*px2-pz2)*vel0*dvdz-v2*px2*dedz+v4*px2*pz2*(dedz-deldz);
	dhdt=v2*(bepos*px2+pz2-4*(epos0-delta0)*v2*px2*pz2);
	printf("dhdpx=%e   dhdpz=%e dhdt=%e\n",dhdpx,dhdpz,dhdt);
	if(cray.rs[cray.nrs-1].x==geo2dv.xmax||cray.rs[cray.nrs-1].x==geo2dv.xmin){
		dtu=-t_x/dhdx;
		pray.rs[pray.nrs].x=cray.rs[cray.nrs-1].x+t_x+dtu*dhdx;
		pray.rs[pray.nrs].z=cray.rs[cray.nrs-1].z+t_z+dtu*dhdz;
		pray.rs[pray.nrs].px= cray.rs[cray.nrs-1].px+t_px+dtu*dhdpx;
		pray.rs[pray.nrs].pz= cray.rs[cray.nrs-1].pz+t_pz+dtu*dhdpz;
		pray.rs[pray.nrs].t = cray.rs[cray.nrs-1].t+dtu*dhdt;
		pray.nrs++;
	}else{
		dtu=-t_z/dhdz;
		pray.rs[pray.nrs].x=cray.rs[cray.nrs-1].x+t_x+dtu*dhdx;
		pray.rs[pray.nrs].z=cray.rs[cray.nrs-1].z+t_z+dtu*dhdz;
		pray.rs[pray.nrs].px= cray.rs[cray.nrs-1].px+t_px+dtu*dhdpx;
		pray.rs[pray.nrs].pz= cray.rs[cray.nrs-1].pz+t_pz+dtu*dhdpz;
		pray.rs[pray.nrs].t = cray.rs[cray.nrs-1].t+dtu*dhdt;
		pray.nrs++;
	}
	printf("dtu = %lf \n",dtu);

//	for(i=0;i<cray.nrs;i++){
//		printf("------------I=%d---------------\n",i);
//		printf("center: X=%lf,Z=%lf,Px=%e,Pz=%e,T=%lf\n",cray.rs[i].x,cray.rs[i].z,cray.rs[i].px,cray.rs[i].pz,cray.rs[i].t);
//		if(i<pray.nrs)
//		  printf("Perturb:X=%lf,Z=%lf,Px=%e,Pz=%e,T=%lf\n",pray.rs[i].x,pray.rs[i].z,pray.rs[i].px,pray.rs[i].pz,pray.rs[i].t);
//	}

	printf("\n");
	printf("=======================================================================\n");
	printf("C Lp: X=%lf,Z=%lf,Px=%1.6e,Pz=%1.6e,T=%lf\n",cray.rs[cray.nrs-1].x,cray.rs[cray.nrs-1].z,cray.rs[cray.nrs-1].px,cray.rs[cray.nrs-1].pz,cray.rs[cray.nrs-1].t);
	printf("P Lp: X=%lf,Z=%lf,Px=%1.6e,Pz=%1.6e,T=%lf\n",pray.rs[pray.nrs-1].x,pray.rs[pray.nrs-1].z,pray.rs[pray.nrs-1].px,pray.rs[pray.nrs-1].pz,pray.rs[pray.nrs-1].t);
	printf("=======================================================================\n");

	FILE *c_path,*p_path;
	c_path=fopen(center_raypath,"w");
	p_path=fopen(perturb_raypath,"w");
	draw_raypath(c_path,cray);
	draw_raypath(p_path,pray);
	fclose(c_path);
	fclose(p_path);
	
	free(vel);
	free(epos);
	free(delta);
	free(dvel);
	free(depos);
	free(ddel);
	free(cray.rs);
	free(pray.rs);
	return 0;
}

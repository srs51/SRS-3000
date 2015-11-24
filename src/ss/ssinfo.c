/*
** ssinfo.c -- DCR 3/21/07
** ========
** Outputs simple information about the contents of ss files.
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <ss.h>
#include <vector.h>
#include <boolean.h>

typedef struct {
  double dTime;
  int nData;
  double dTotMass,*dMass,*dRadius;
  VECTOR *vPos,*vVel,*vSpin;
  int *iColor,*iOrgIdx;
  } AllData;

typedef struct {
  double dMin,dMax;
  int iMinIdx,iMaxIdx,iMinOrg,iMaxOrg;
  } sLimits;

typedef struct {
  VECTOR vMin,vMax,vCom;
  sLimits mag,mag_com;
  } vLimits;

typedef struct {
  int iMin,iMax,iMinIdx,iMaxIdx,iMinOrg,iMaxOrg;
  } iLimits;

static int read_ssfile(const char *achFile,AllData *a)
{
	SSIO ssio;
	SSHEAD h;
	SSDATA d;
	int i;

	if (ssioOpen(achFile,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open ss file for reading.\n");
		return 1;
		}

	if (ssioHead(&ssio,&h) || h.n_data <= 0) {
		(void) fprintf(stderr,"Corrupt header.\n");
		(void) ssioClose(&ssio);
		return 1;
		}

	switch(h.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"Reduced ss format not supported.\n");
		ssioClose(&ssio);
		return 1;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio);
		return 1;
		}

	a->dTime = h.time;
	a->nData = h.n_data;

	a->dMass = (double *) realloc(a->dMass,a->nData*sizeof(double));
	assert(a->dMass);
	a->dRadius = (double *) realloc(a->dRadius,a->nData*sizeof(double));
	assert(a->dRadius);
	a->vPos = (VECTOR *) realloc(a->vPos,a->nData*sizeof(VECTOR));
	assert(a->vPos);
	a->vVel = (VECTOR *) realloc(a->vVel,a->nData*sizeof(VECTOR));
	assert(a->vVel);
	a->vSpin = (VECTOR *) realloc(a->vSpin,a->nData*sizeof(VECTOR));
	assert(a->vSpin);
	a->iColor = (int *) realloc(a->iColor,a->nData*sizeof(int));
	assert(a->iColor);
	a->iOrgIdx = (int *) realloc(a->iOrgIdx,a->nData*sizeof(int));
	assert(a->iOrgIdx);

	a->dTotMass = 0.0;

	for (i=0;i<a->nData;i++) {
		if (ssioData(&ssio,&d)) {
			(void) fprintf(stderr,"Corrupt data (particle %i).\n",i);
			(void) ssioClose(&ssio);
			return 1;
			}
		a->dTotMass += (a->dMass[i] = d.mass);
		a->dRadius[i] = d.radius;
		COPY_VEC(d.pos,a->vPos[i]);
		COPY_VEC(d.vel,a->vVel[i]);
		COPY_VEC(d.spin,a->vSpin[i]);
		a->iColor[i] = d.color;
		a->iOrgIdx[i] = d.org_idx;
		}

	return 0;
	}

static void get_scalar_limits(const AllData *a,const double d[],sLimits *lim)
{
	int i;

	lim->dMin = HUGE_VAL;
	lim->dMax = -HUGE_VAL;
	lim->iMinIdx = lim->iMaxIdx = INT_MAX;
	lim->iMinOrg = lim->iMaxOrg = INT_MAX;

	for (i=0;i<a->nData;i++) {
		if (d[i] < lim->dMin) {
			lim->dMin = d[i];
			lim->iMinIdx = i;
			lim->iMinOrg = a->iOrgIdx[i];
			}
		if (d[i] > lim->dMax) {
			lim->dMax = d[i];
			lim->iMaxIdx = i;
			lim->iMaxOrg = a->iOrgIdx[i];
			}
		}
	}

static void get_vector_limits(const AllData *a,VECTOR v[],vLimits *lim,BOOLEAN bDoCom)
{
	VECTOR vTmp;
	double *dMag;
	int i;

	SET_VEC(lim->vMin,HUGE_VAL,HUGE_VAL,HUGE_VAL);
	SET_VEC(lim->vMax,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL);
	if (bDoCom) {
		ZERO_VEC(lim->vCom);
		}

	for (i=0;i<a->nData;i++) {
		if (v[i][X] < lim->vMin[X])
			lim->vMin[X] = v[i][X];
		if (v[i][Y] < lim->vMin[Y])
			lim->vMin[Y] = v[i][Y];
		if (v[i][Z] < lim->vMin[Z])
			lim->vMin[Z] = v[i][Z];
		if (v[i][X] > lim->vMax[X])
			lim->vMax[X] = v[i][X];
		if (v[i][Y] > lim->vMax[Y])
			lim->vMax[Y] = v[i][Y];
		if (v[i][Z] > lim->vMax[Z])
			lim->vMax[Z] = v[i][Z];
		if (bDoCom) {
			COPY_VEC(v[i],vTmp);
			SCALE_VEC(vTmp,a->dMass[i]);
			ADD_VEC(lim->vCom,vTmp,lim->vCom);
			}
		}

	if (bDoCom) {
		if (a->dTotMass > 0.0) {
			NORM_VEC(lim->vCom,a->dTotMass);
			}
		else {
			NORM_VEC(lim->vCom,a->nData);
			}
		}

	dMag = (double *) malloc(a->nData*sizeof(double));
	assert(dMag != NULL);

	for (i=0;i<a->nData;i++) {
		dMag[i] = MAG(v[i]);
		}

	get_scalar_limits(a,dMag,&lim->mag);

	if (bDoCom) {
		for (i=0;i<a->nData;i++) {
			SUB_VEC(v[i],lim->vCom,vTmp);
			dMag[i] = MAG(vTmp);
			}
		get_scalar_limits(a,dMag,&lim->mag_com);
		}

	free((void *) dMag);
	}

static void get_int_limits(const AllData *a,const int iVal[],iLimits *lim)
{
	int i;

	lim->iMin = INT_MAX;
	lim->iMax = INT_MIN;
	lim->iMinIdx = lim->iMaxIdx = INT_MAX;
	lim->iMinOrg = lim->iMaxOrg = INT_MAX;

	for (i=0;i<a->nData;i++) {
		if (iVal[i] < lim->iMin) {
			lim->iMin = iVal[i];
			lim->iMinIdx = i;
			lim->iMinOrg = a->iOrgIdx[i];
			}
		if (iVal[i] > lim->iMax) {
			lim->iMax = iVal[i];
			lim->iMaxIdx = i;
			lim->iMaxOrg = a->iOrgIdx[i];
			}
		}
	}

static void show_scalar_limits(const sLimits *lim,const char *achType,double dScale,const char *achUnits)
{
	(void) printf("%s range: %g (%g %s) [%i,%i] to %g (%g %s) [%i,%i]\n",achType,
				  lim->dMin,lim->dMin*dScale,achUnits,lim->iMinIdx,lim->iMinOrg,
				  lim->dMax,lim->dMax*dScale,achUnits,lim->iMaxIdx,lim->iMaxOrg);
	}

static void show_vector_limits(const vLimits *lim,const char *achType,double dScale,const char *achUnits,BOOLEAN bDoCom)
{
	(void) printf("%s X range: %g (%g %s) to %g (%g %s)\n",achType,
				  lim->vMin[X],lim->vMin[X]*dScale,achUnits,
				  lim->vMax[X],lim->vMax[X]*dScale,achUnits);
	(void) printf("%s Y range: %g (%g %s) to %g (%g %s)\n",achType,
				  lim->vMin[Y],lim->vMin[Y]*dScale,achUnits,
				  lim->vMax[Y],lim->vMax[Y]*dScale,achUnits);
	(void) printf("%s Z range: %g (%g %s) to %g (%g %s)\n",achType,
				  lim->vMin[Z],lim->vMin[Z]*dScale,achUnits,
				  lim->vMax[Z],lim->vMax[Z]*dScale,achUnits);

	if (bDoCom)
		(void) printf("%s com: %g,%g,%g (%g,%g,%g %s)\n",achType,
					  lim->vCom[X],lim->vCom[Y],lim->vCom[Z],
					  lim->vCom[X]*dScale,lim->vCom[Y]*dScale,lim->vCom[Z]*dScale,achUnits);
	}

static void show_int_limits(const iLimits *lim,const char *achType)
{
	(void) printf("%s range: %i [%i,%i] to %i [%i,%i]\n",achType,
				  lim->iMin,lim->iMinIdx,lim->iMinOrg,
				  lim->iMax,lim->iMaxIdx,lim->iMaxOrg);
	}				  

int main(int argc,char *argv[])
{
	AllData all;
	sLimits sLim;
	vLimits vLim;
	iLimits iLim;
	int i;

	if (argc < 2) {
		(void) fprintf(stderr,"Usage: %s ss-file [ ss-file ... ]\n",argv[0]);
		return 1;
		}

	all.dMass = all.dRadius = NULL;
	all.vPos = all.vVel = all.vSpin = NULL;
	all.iColor = all.iOrgIdx = NULL;

	for (i=1;i<argc;i++) {
		(void) printf("%s:\n",argv[i]);
		if (read_ssfile(argv[i],&all) == 0) {
			(void) printf("Time = %g (%g yr)\n",all.dTime,all.dTime*T_SCALE/SID_YR);
			(void) printf("Number of particles = %i\n",all.nData);
			(void) printf("Total mass = %g (%g kg)\n",all.dTotMass,all.dTotMass*M_SCALE);
			get_scalar_limits(&all,all.dMass,&sLim);
			show_scalar_limits(&sLim,"Mass",M_SCALE,"kg");
			get_scalar_limits(&all,all.dRadius,&sLim);
			show_scalar_limits(&sLim,"Radius",0.001*L_SCALE,"km");
			get_vector_limits(&all,all.vPos,&vLim,TRUE);
			show_vector_limits(&vLim,"Position",0.001*L_SCALE,"km",TRUE);
			show_scalar_limits(&vLim.mag,"Pos abs mag",0.001*L_SCALE,"km");
			show_scalar_limits(&vLim.mag_com,"Pos com mag",0.001*L_SCALE,"km");
			get_vector_limits(&all,all.vVel,&vLim,TRUE);
			show_vector_limits(&vLim,"Velocity",V_SCALE,"m/s",TRUE);
			show_scalar_limits(&vLim.mag,"Vel abs mag",V_SCALE,"m/s");
			show_scalar_limits(&vLim.mag_com,"Vel com mag",V_SCALE,"m/s");
			get_vector_limits(&all,all.vSpin,&vLim,FALSE);
			show_vector_limits(&vLim,"Spin",1.0/T_SCALE,"rad/s",FALSE);
			show_scalar_limits(&vLim.mag,"Spin mag",1.0/T_SCALE,"rad/s");
			get_int_limits(&all,all.iColor,&iLim);
			show_int_limits(&iLim,"Color");
			get_int_limits(&all,all.iOrgIdx,&iLim);
			show_int_limits(&iLim,"Original index");
			}
		(void) printf("\n");
		}

	if (all.iOrgIdx == NULL)
		return 1;

	free((void *) all.iOrgIdx);
	free((void *) all.iColor);
	free((void *) all.vSpin);
	free((void *) all.vVel);
	free((void *) all.vPos);
	free((void *) all.dRadius);
	free((void *) all.dMass);

	return 0;
	}

/* ssinfo.c */

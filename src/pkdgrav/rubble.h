#ifndef _RUBBLE_H
#define _RUBBLE_H

#ifdef RUBBLE_ZML

#include "pkd.h"
#include "ssdefs.h" /* in turn includes ssio.h */

#define M_SUN	1.9891e30		/* One solar mass in kg, from ss_core's ss.h */
#define AU		1.49597870e11	/* One AU in metres, from ss_core's ss.h */

#define TWO_PI (2.0*M_PI)

double SQ(double);
double CUBE(double);
int SGN(double);

#define SQ(x)    ((x)*(x))     				/* Square of x */
#define CUBE(x)  ((x)*(x)*(x))			       	/* Cube of x */
#define SGN(x)   ((x) < 0.0 ? (-1) : (x) > 0.0 ? 1 : 0)	/* Signum of x */

#define MAX_NPART 100 /* number of particles for target (projectile may have fewer) */
#define SUPER_NPART 2500 /*DEBUG 07.17.04*/
/* DEBUG: not all values might not work test before using 04.19.04 */
#define RUB_MAX_NUM_EVENTS 50 /* temporary? */
#define RUB_BASE_COLOR 4 /* colors 4 to (4 + RUB_MAX_NUM_EVENTS - 1) reserved for rubble piles */

#define RUB_FORCED_NONE 0
#define RUB_FORCED_BOUNCE 1
#define RUB_FORCED_MERGE 2

typedef struct {
	int iColor;
	double dTStartMergePhase;
	double dTEndRubblePhase;
	} RUB_CLOCK;

struct rubEvents {
	int nEvents;
	RUB_CLOCK rc[RUB_MAX_NUM_EVENTS];
	};

/* dustbin parameters */

#define DUST_BINS_MAX_NUM 40

/* following stored in msr->params */

typedef struct { /* see ss.par for legend */
	int nDustBins;
	int iDustBinsApplyInt;
	int iRubNumDynToBounce;
	int iRubNumDynToMerge;
  int iDustBinsVelDispOpt;
	double dDustBinsInner;
	double dDustBinsOuter;
	double dDustBinsScaleHeight;
	double dDustBinsInitSigma;
	double dDustBinsInitAlpha;
	double dDustBinsWidth; /* computed */
	double dRubMinMass;
	} DUST_BINS_PARAMS;

/* following stored in msr */

typedef struct { /* there will be nDustBins of these */
	double dMass;
	double dVolume;
#ifdef ORIGIN_HISTOGRAM
	FLOAT origin_bins[NUM_ORIGIN_BINS];
#endif /* ORIGIN_HISTOGRAM */ 
	} DustBins;

/* pkd routine prototypes */
double pkdDustBinsInclMax(PKD pkd);
double pkdDustBinsInclAvg(PKD pkd);
void pkdDustBinsGetMass(PKD pkd,DUST_BINS_PARAMS *DB,DustBins pDustBins[],
                        double dt,double M,double pDustBinsMassLoss[]);
void pkdDustBinsApply(PKD pkd,double M,double pMassIncrFrac[],int nBins,
		      double dt,DustBins pDustBins[],double dVdisp
#ifdef ORIGIN_HISTOGRAM
					  ,FLOAT aOriginBins[][NUM_ORIGIN_BINS]
#endif /* ORIGIN_HISTOGRAM */
					  );
void pkdRubbleResetColFlag(PKD pkd);
int pkdRubbleCheckForKDKRestart(PKD pkd);
void pkdRubbleStep(PKD pkd,double dMaxStep,double dMinStep);
void pkdRubCleanup(PKD pkd,int iColor,DUST_BINS_PARAMS *DB,double M,
                   DustBins aDustBins[],DustBins *pDustBinsTrash);
void pkdRubInterpCleanup(PKD pkd,DUST_BINS_PARAMS *DB,double M,int iOrder,
                         int *iBin,DustBins *DustBinsInterp);

/* from vector.h -- DCR 94-04-18 */
/* zml added to rubble.h 03-04-01 */

/*** Three dimensions... ***/

#ifndef N_DIM
#  define N_DIM 3
#endif

/*** Vector and matrix type definitions ***/

typedef double VECTOR[N_DIM];
typedef double MATRIX[N_DIM][N_DIM];

/*** Vector component definitions ***/

#define X 0
#define Y 1
#define Z 2

/*** Macro definitions ***/

/* Assigns a value (x,y,z) to vector v */

#define SET_VEC(v, x, y, z) {\
								 (v)[X] = (x);\
								 (v)[Y] = (y);\
								 (v)[Z] = (z);\
							 }

/* Assigns zero to vector v */

#define ZERO_VEC(v) SET_VEC((v), 0, 0, 0)

/* Copies vector v1 to vector v2 */

#define COPY_VEC(v1, v2) {\
							  (v2)[X] = (v1)[X];\
							  (v2)[Y] = (v1)[Y];\
							  (v2)[Z] = (v1)[Z];\
						  }

/* Adds vectors v1 & v2 and puts the result in vector v */

#define ADD_VEC(v1, v2, v) {\
								(v)[X] = (v1)[X] + (v2)[X];\
								(v)[Y] = (v1)[Y] + (v2)[Y];\
								(v)[Z] = (v1)[Z] + (v2)[Z];\
							}

/* Subtracts vector v2 from vector v1 and puts the result in vector v */

#define SUB_VEC(v1, v2, v) {\
								(v)[X] = (v1)[X] - (v2)[X];\
								(v)[Y] = (v1)[Y] - (v2)[Y];\
								(v)[Z] = (v1)[Z] - (v2)[Z];\
							}

/* Multiplies vector v by scalar a */

#define SCALE_VEC(v, a) {\
							 (v)[X] *= (a);\
							 (v)[Y] *= (a);\
							 (v)[Z] *= (a);\
						 }

/* Divides vector v by scalar a */

#define NORM_VEC(v, a) {\
							double _n = 1.0 / (a);\
							SCALE_VEC((v), _n);\
						}

/* Returns dot product of vectors v1 & v2 */

#define DOT(v1, v2) ((v1)[X] * (v2)[X] + (v1)[Y] * (v2)[Y] + (v1)[Z] * (v2)[Z])

/* Returns cross product of vectors v1 & v2 in vector v */

#define CROSS(v1, v2, v) {\
							  (v)[X] = (v1)[Y] * (v2)[Z] - (v1)[Z] * (v2)[Y];\
							  (v)[Y] = (v1)[Z] * (v2)[X] - (v1)[X] * (v2)[Z];\
							  (v)[Z] = (v1)[X] * (v2)[Y] - (v1)[Y] * (v2)[X];\
						  }

/* Returns square magnitude of vector v */

#define MAG_SQ(v) (DOT((v), (v)))

/* Returns magnitude of vector v */

#define MAG(v) (sqrt(MAG_SQ(v)))

/* Zeroes 3x3 matrix a */

#define ZERO_MAT(a) {\
						 ZERO_VEC((a)[X]);\
						 ZERO_VEC((a)[Y]);\
						 ZERO_VEC((a)[Z]);\
					 }

/* Makes matrix a the unit matrix */

#define UNIT_MAT(a) {\
						 ZERO_MAT(a);\
						 (a)[X][X] = (a)[Y][Y] = (a)[Z][Z] = 1;\
					 }

/* Copies matrix a to matrix b */

#define COPY_MAT(a, b) {\
							COPY_VEC((a)[X], (b)[X]);\
							COPY_VEC((a)[Y], (b)[Y]);\
							COPY_VEC((a)[Z], (b)[Z]);\
						}

/* Multiplies matrix a by scalar b */

#define SCALE_MAT(a, b) {\
							 SCALE_VEC((a)[X], (b));\
							 SCALE_VEC((a)[Y], (b));\
							 SCALE_VEC((a)[Z], (b));\
						 }

/* from boolean header file dcr 98-09-16 */
/* zml added to rubble.h 03-04-01 */

#ifndef BOOLEAN
#define BOOLEAN int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

/* from rpu header file dcr 98-09-23*/
/* zml added to rubble.h 03-04-01 */

#define MAJOR(prp) ((prp)->axis_ord[0])
#define INTER(prp) ((prp)->axis_ord[1])
#define MINOR(prp) ((prp)->axis_ord[2])

typedef struct {
	/* traditional quantities */
	double mass;
	double radius;
	VECTOR pos;
	VECTOR vel;
	VECTOR spin;
	/* special quantities */
	int agg_id; /*DEBUG not used*/
	MATRIX inertia;
	MATRIX axes;
	VECTOR ang_mom; /*DEBUG not used*/
	VECTOR moments;
	VECTOR axis_len; /* semi-axes */
	int axis_ord[N_DIM];
	double density;
	double packing;
	/* particle data */
	int n_particles;
	SSDATA *data;
	} RUBBLE_PILE;

#endif /* RUBBLE_ZML */

#endif /* _RUBBLE_H */
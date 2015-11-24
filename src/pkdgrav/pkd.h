#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <sys/time.h>
#include <sys/resource.h>
#include "mdl.h"
#include "floattype.h"

#ifdef GASOLINE
#include "cooling.h"
#endif
#include "rotbar.h"

#ifdef SPRINGS
#include "springs.h"
#endif

#ifdef DEM
#include "dem.h"
#endif

#ifdef SLIDING_PATCH
#include "patch.h"
#endif

#ifdef JOHNNY
#include "delaunay.h"
#include "helio.h"
#endif

/*
 ** The following sort of definition should really be in a global
 ** configuration header file -- someday...
 */

#if defined(GASOLINE) || defined(ROT_FRAME) || defined(SIMPLE_GAS_DRAG) || defined(SPRINGS) || defined(DEM)
#define NEED_VPRED
#endif

/* this too... */

/* (note bVWarnings still applies) */

#define INTERNAL_WARNINGS 1000000 /* 0=none,1=always,>1=every n-th occurence (plus first) */

#define CID_TOP			0
#define CID_PARTICLE	0
#define CID_CELL		1

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define MAX_TIMERS		10

typedef struct particle {
	int iOrder;
	unsigned int iActive;  
	int iRung;
	int cpStart;
	FLOAT fWeight;
	FLOAT fMass;
	FLOAT fSoft;
#ifdef CHANGESOFT
	FLOAT fSoft0;
#endif
	FLOAT r[3];
	FLOAT v[3];
	FLOAT a[3];
	FLOAT fPot;
	FLOAT fBall2;
	FLOAT fDensity;
	FLOAT dt;			/* a time step suggestion */
	FLOAT dtGrav;		/* suggested 1/dt^2 from gravity */
#ifdef SLIDING_PATCH
	FLOAT dPy;		/* Canonical momentum for Hill eqn. */
#endif
#ifdef SUPERCOOL
	FLOAT vMean[3];
#endif
#ifdef COLORCODE
	FLOAT fColor;
#endif
	FLOAT fBallMax;		/* SPH 2h Max value (now also used for COLLISIONS searches!) */
#ifdef GASOLINE
	FLOAT uPred;		/* predicted thermal energy */
	FLOAT PoverRho2;	/* P/rho^2 */
	FLOAT u;	        /* thermal energy */ 
	FLOAT c;			/* sound speed */
	FLOAT mumax;		/* sound speed like viscosity term */
	FLOAT PdV;	        /* P dV heating (includes shocking) */
#ifdef PDVDEBUG
	FLOAT PdVvisc;		/* P dV from shock (testing) */
	FLOAT PdVpres;		/* P dV from adiabatic compression (testing) */
#endif
	FLOAT divv;		/* Balsara viscosity reduction */
	FLOAT curlv[3];         
	FLOAT BalsaraSwitch;
#ifdef SHOCKTRACK
        FLOAT aPres[3];
        FLOAT ShockTracker;     /* Shock tracker */
        FLOAT divrhov;          /* debug */
        FLOAT gradrho[3];       /* debug */
#endif
/*	FLOAT fDensSave;*/	/* Used by diagnostic DensCheck funcs */
#ifndef NOCOOLING
	FLOAT uDot;			/* Rate of change of u -- for predicting u */
	COOLPARTICLE CoolParticle;  /* Abundances and any other cooling internal variables */
#endif
	FLOAT fMetals;
	FLOAT fTimeForm;
#ifdef SIMPLESF
    FLOAT fMassStar;
	FLOAT fESN;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	int iGasOrder;		/* gas from which star formed */
#endif
#ifdef STARFORM
	FLOAT fESNrate;
	FLOAT fMSN;
	FLOAT fNSN;           
	FLOAT fMOxygenOut;
	FLOAT fMIronOut;
	FLOAT fMFracOxygen;
	FLOAT fMFracIron;
	FLOAT fSNMetals;
	FLOAT fNSNtot;
        FLOAT fTimeCoolIsOffUntil;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	FLOAT fMassForm;	/* record original mass of star */
	int iGasOrder;		/* gas from which star formed */
#endif
#endif  /* GASOLINE */
#ifdef COLLISIONS
	int iOrgIdx;		/* for tracking of mergers, aggregates etc. */
	FLOAT w[3];			/* spin vector */
	int iColor;			/* handy color tag */
	double dtCol;		/* time to next encounter or collision */
	int iOrderCol;		/* neighbour or collider iOrder */
	double dtPrevCol;	/* time of previous collision */
	int iPrevCol;		/* iOrder of previous collider */
	int bTinyStep;		/* flag for imminent collapse */
        FLOAT fOverlap; /* max fractional overlap with other particle */
#endif /* COLLISIONS */
#ifdef SPECIAL_PARTICLES
	int bGhostExclude;	/* particle not included in ghost cells */
#endif /* SPECIAL_PARTICLES */ 
#ifdef SLIDING_PATCH
        int bAzWrap;        /* flag set on azimuthal boundary wrap */
#endif
#ifdef SPRINGS
  SPRING springs[MAX_NUM_SPRINGS_PER_PARTICLE];
#endif
#ifdef DEM
  FLOAT wDot[3]; /* rate of change of angular velocity */
  FLOAT wPred[3]; /* predicted angular velocity (time centered) */
  FLOAT dPressure; /* sum of all normal forces on particle */
  DEM_ELEMENT overlaps[MAX_NUM_OVERLAPS_PER_PARTICLE];
#ifdef WALLS
  DEM_ELEMENT walloverlaps[MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS];
#ifdef WALLS_REACT
  double dZForceOnWalls;
#endif /* WALLS_REACT */
#endif /* WALLS */
  DEM_STATS_PARTICLE DEMStats;
#endif /* DEM */
#if defined(SPRINGS) || defined(DEM)
  FLOAT dDeltaAccel[3]; /* need combiner cache? */
#endif
#ifdef CHARGE
  double dCharge; /* in units of electron charge (1.602e-19 C) */
#endif /* CHARGE */
#ifdef NEED_VPRED
  FLOAT vPred[3];		/* predicted velocity (time centered) */
#endif
#ifdef AGGS
  /*
  ** Position of particle in principal frame of the aggregate
  ** (normally).  We temporarily store the COM frame position here
  ** during the process of computing the aggregate parameters.
  */
  FLOAT r_body[3];
#endif
#if defined(AGGS) || defined(WALLS)
  FLOAT omega_v[3]; /* 2nd-order expansion for collision prediction */
#endif
#if defined(AGGS) && defined(DEM)
  FLOAT fMassAgg; /* mass of the aggregate to which the particle belongs */
#endif
#if defined (RUBBLE_ZML) || defined (COLLMOD_ZML)
  double dDustMass;	// predicted mass increase from dust
  int iBin;			// dust bin that planetesimal is in
  int bMayCollide;	// true if planetesimal is predicted to collide in this step
#endif
#ifdef COLLMOD_ZML
  int ctimer;		// time assigned for suppressing debris collisions
#endif
#ifdef ORIGIN_HISTOGRAM
	FLOAT origin_bins[NUM_ORIGIN_BINS];
#endif /* ORIGIN_HISTOGRAM */
#ifdef AGGS_IN_PATCH
  /*
  ** For this special case, we need to keep track of particle
  ** quantities with and without boundary wrapping.  The traditional r
  ** & v vectors will contain the *wrapped* quantities (so they are
  ** defined within the patch boundaries) and the following will
  ** contain the unwrapped quantities (so we can more easily compute
  ** things like the aggregate center of mass location, since the
  ** particles will be close to the aggregate center, but possibly
  ** outside the patch boundaries).
  */
  FLOAT r_unwrapped[3]; /* note: z component same as r[2] */
  FLOAT v_unwrapped[3]; /* note: x & z same as v[0] & v[2] */
  FLOAT v_patch[3]; /* needed for collision prediction & backdrifting */
#endif /* AGGS_IN_PATCH */
#ifdef GR_DRAG
  int bNoKickNoDrift;
  double dLastPericenter;
  double dEntryR2;
#endif
#ifdef JOHNNY
  int bInKepler;
  int bRepeat;
  struct delaunay sOrbElem;
  double dMotionPerStep;
#endif
} PARTICLE;

/* Active Type Masks */

/* Active: -- eg. Calculate new acceleration, PdV, etc... for this particle */
#define TYPE_ACTIVE            (1<<0)
/* In the Tree: */
#define TYPE_TREEACTIVE        (1<<1)
/* Gather to/Scatter from this particle with in smooths: */
#define TYPE_SMOOTHACTIVE      (1<<2)
/* Smooth has processed this particle */
#define TYPE_SMOOTHDONE        (1<<3)

/* Types used for Fast Density only (so far) */
/* Sum Fast Density on this particle */
#define TYPE_DensACTIVE        (1<<4)
/* Neighbour of ACTIVE (incl. ACTIVE): */
#define TYPE_NbrOfACTIVE       (1<<5)
/* Potential Scatter Neighbour */
#define TYPE_Scatter           (1<<6)
/* Density set to zero already */
#define TYPE_DensZeroed        (1<<7)

/* Particle Type Masks */

#define TYPE_GAS               (1<<8)
#define TYPE_DARK              (1<<9)
#define TYPE_STAR              (1<<10)
#define TYPE_SUPERCOOL         (1<<11)

/* Particle marked for deletion.  Will be deleted in next
   msrAddDelParticles(); */
#define TYPE_DELETED           (1<<12)

#define TYPE_PHOTOGENIC        (1<<13)
#define TYPE_SINK              (1<<14)

/* Combination Masks */
#define TYPE_ALLACTIVE			(TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE)
#define TYPE_ALL				(TYPE_GAS|TYPE_DARK|TYPE_STAR)

/* Type Macros */
int TYPEQueryACTIVE      ( PARTICLE *a );
int TYPEQueryTREEACTIVE  ( PARTICLE *a );
int TYPEQuerySMOOTHACTIVE( PARTICLE *a );
int TYPETest  ( PARTICLE *a, unsigned int mask );
int TYPEFilter( PARTICLE *a, unsigned int filter, unsigned int mask );
int TYPESet   ( PARTICLE *a, unsigned int mask );
int TYPEReset ( PARTICLE *a, unsigned int mask );
/* This retains Particle Type and clears all flags: */
int TYPEClearACTIVE( PARTICLE *a ); 

/* Warning: This erases Particle Type */
int TYPEClear( PARTICLE *a ); 

#define TYPEQueryACTIVE(a)       ((a)->iActive & TYPE_ACTIVE)
#define TYPEQueryTREEACTIVE(a)   ((a)->iActive & TYPE_TREEACTIVE)
#define TYPEQuerySMOOTHACTIVE(a) ((a)->iActive & TYPE_SMOOTHACTIVE)
#define TYPETest(a,b)            ((a)->iActive & (b))
#define TYPEFilter(a,b,c)        (((a)->iActive & (b))==(c))
#define TYPESet(a,b)             ((a)->iActive |= (b))
#define TYPEReset(a,b)           ((a)->iActive &= (~(b)))
#define TYPEClearACTIVE(a)       ((a)->iActive &= (TYPE_ALL|TYPE_SUPERCOOL))
#define TYPEClear(a)             ((a)->iActive = 0)

#ifdef RUBBLE_ZML
/* RUBBLE_ZML (and ORIGIN_HISTOGRAM) puts extra stuff in the checkpoint file */
	#ifdef ORIGIN_HISTOGRAM
		#define CHECKPOINT_VERSION 82
	#else
		#define CHECKPOINT_VERSION 81
	#endif
#elif COLLMOD_ZML
	#ifdef ORIGIN_HISTOGRAM
		#define CHECKPOINT_VERSION 92
	#else
		#define CHECKPOINT_VERSION 91
	#endif
#else
	#define CHECKPOINT_VERSION 8
#endif /* RUBBLE_ZML, COLLMOD_ZML */

/*
** NOTE: the #ifdef's below, in the checkpoint particle structure,
** imply that any attempt to read a checkpoint that was written with
** different compile flags will likely fail spectacularly.  You have
** been warned!
*/

typedef struct chkParticle {
	int iOrder;
	FLOAT fMass;
	FLOAT fSoft;
	FLOAT r[3];
	FLOAT v[3];
#ifdef GASOLINE
	FLOAT u;
	FLOAT fMetals;
#ifndef NOCOOLING
	COOLPARTICLE CoolParticle;
#endif
#ifdef STARFORM
        FLOAT fTimeCoolIsOffUntil;
	FLOAT fTimeForm;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	FLOAT fMassForm;	/* record original mass of star */
	FLOAT fDenForm;
	FLOAT fNSN;
	FLOAT fMFracOxygen;
	FLOAT fMFracIron;
        int iGasOrder;
#endif
#ifdef SIMPLESF
	FLOAT fMassStar;
	FLOAT fTimeForm;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	FLOAT fDenForm;
	int iGasOrder;
#endif
#endif
#ifdef COLLISIONS
	int iOrgIdx; /* added for version 7 */
	FLOAT w[3];
	int iColor;
#ifdef SPRINGS
  SPRING springs[MAX_NUM_SPRINGS_PER_PARTICLE];
#endif /* SPRINGS */
#ifdef DEM
  FLOAT wDot[3]; /* rate of change of angular velocity */
  FLOAT wPred[3]; /* predicted angular velocity (time centered) */
  DEM_ELEMENT overlaps[MAX_NUM_OVERLAPS_PER_PARTICLE];
#ifdef WALLS
  DEM_ELEMENT walloverlaps[MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS];
#endif /* WALLS */
#endif /* DEM */
#endif /* COLLISIONS */
#ifdef ORIGIN_HISTOGRAM
  FLOAT origin_bins[NUM_ORIGIN_BINS];
#endif /* ORIGIN_HISTOGRAM */ 
#ifdef GR_DRAG
  int bNoKickNoDrift;
  double dLastPericenter;
  double dEntryR2;
#endif
#ifdef JOHNNY
  int bInKepler;
#endif /* JOHNNY */
  } CHKPART;

typedef struct bndBound {
  FLOAT fMin[3];
  FLOAT fMax[3];
  } BND;

struct pkdCalcCellStruct {
	double Qxx,Qyy,Qzz,Qxy,Qxz,Qyz;
	/*
	 ** Reduced multipole moments for l>2 !!!
	 */
	double Oxxx,Oxyy,Oxxy,Oyyy,Oxxz,Oyyz,Oxyz;
	double Oxzz, Oyzz, Ozzz;
	double Hxxxx,Hxyyy,Hxxxy,Hyyyy,Hxxxz,Hyyyz,Hxxyy,Hxxyz,Hxyyz;
	double Hxxzz, Hxyzz, Hxzzz, Hyyzz, Hyzzz, Hzzzz;
	double Bmax,B2,B3,B4,B5,B6;
	};

typedef struct kdNode {
	int iDim;
	double fSplit;
	BND bnd;
	BND bndBall;	/* Bound including fBall*(1+changemax) */
	int pLower;		/* also doubles as thread id for the LTT */
	int pUpper;		/* pUpper < 0 indicates no particles in tree! */
	int iLower;
	int iUpper;
	double fMass;
	double fSoft;
	FLOAT r[3];
	struct pkdCalcCellStruct mom;
	double fOpen2;
	} KDN;

typedef struct ilPart {
	double m,h;
	double x,y,z;
	} ILP;

typedef struct ilCellSoft {
	double m,h;
	double x,y,z;
	double xx,yy,zz,xy,xz,yz;
	} ILCS;

/*
 ** moment tensor components.
 */
typedef struct ilCellNewt {
	double m;
	double x,y,z;
	double xx,yy,xy,xz,yz;
	double zz;
	double xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	double xzz,yzz,zzz;
	double xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	double xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} ILCN;

/* IBM brain damage */
#undef hz

typedef struct ewaldTable {
	double hx,hy,hz;
	double hCfac,hSfac;
	} EWT;

typedef struct pkdContext {
	MDL mdl;
	int idSelf;
	int nThreads;
	int nStore;
	int nRejects;
	int nLocal;
	int nActive;
	int nTreeActive;
	int nSmoothActive;
	int nDark;
	int nGas;
	int nStar;
	int nMaxOrderDark;
	int nMaxOrderGas;
	int nBucket;
	int nLevels;
	int nSplit;
	int nNodes;
	int iExtraBucket;
	int iOrder;
	int iFreeCell;
	int iRoot;
	FLOAT fPeriod[3];
	int *piLeaf;
	KDN *kdTop;
	KDN *kdNodes;
	PARTICLE *pStore;
        double duTFac;
	/*
	 ** gravitational interaction lists
	 */
	int nMaxPart;
	int nMaxCellSoft;
	int nMaxCellNewt;
	int nPart;
	int nCellSoft;
	int nCellNewt;
	int nSqrtTmp;
	ILP *ilp;
	ILCS *ilcs;
	ILCN *ilcn;
	double *sqrttmp;
	double *d2a;
	/*
	 ** Ewald summation setup.
	 */
	ILCN ilcnRoot;
	int nMaxEwhLoop;
	int nEwhLoop;
	EWT *ewt;
	/*
	 ** Timers stuff.
	 */
	struct timer {
		double sec;
		double stamp;
		double system_sec;
		double system_stamp;
		double wallclock_sec;
		double wallclock_stamp;
		int iActive;
		} ti[MAX_TIMERS];
#ifdef GASOLINE
	/* 
	 ** Cooling 
	 */
#ifndef NOCOOLING
	COOL *Cool;
#endif
#endif
#ifdef COLLISIONS
  int bRepel; /* reverse gravity enabled? */
  FLOAT dRepelFac; /* always positive (As of July '09) */
#ifdef SPRINGS
  int bReadSpringsData;
#endif
#ifdef DEM
  int bReadDEMData;
#endif
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
  /* extra quantities needed when dealing with sliding patches... */
  double dTime;
  PATCH_PARAMS *PP;
#endif
        ROTBAR  rotbar;
	} * PKD;

int pkdIsGas(PKD,PARTICLE *);
#define pkdIsGas( pkd, pTMP) TYPETest( (pTMP), TYPE_GAS )
int pkdIsDark(PKD,PARTICLE *);
#define pkdIsDark( pkd, pTMP) TYPETest( (pTMP), TYPE_DARK )
int pkdIsStar(PKD,PARTICLE *);
#define pkdIsStar( pkd, pTMP) TYPETest( (pTMP), TYPE_STAR )

int pkdIsGasByOrder(PKD pkd,PARTICLE *p);
int pkdIsDarkByOrder(PKD pkd,PARTICLE *p);
int pkdIsStarByOrder(PKD pkd,PARTICLE *p);


typedef struct CacheStatistics {
	double dpNumAccess;
	double dpMissRatio;
	double dpCollRatio;
	double dpMinRatio;
	double dcNumAccess;
	double dcMissRatio;
	double dcCollRatio;
	double dcMinRatio;
	} CASTAT;

/* JW: */

#define GASMODEL_UNSET -1
enum GasModel {
	GASMODEL_ADIABATIC, 
	GASMODEL_ISOTHERMAL, 
	GASMODEL_COOLING, 
	GASMODEL_GLASS
	}; 

#define PKD_ORDERTEMP	256

#define pkdRoot(iCell,id)\
{\
	id = -1;\
	iCell = ROOT;\
	}

#define pkdIsRoot(iCell,id)		((id==-1)?((iCell==ROOT)?1:0):0)

/*
 * There is now a slight inefficency here when going from the top tree to
 * a node tree in that we visit the root cell twice (once as the leaf of the
 * top tree and once as the root of the node tree).  This is necessary to
 * check if the root cell is a bucket.
 */

#define pkdLower(pkd,iCell,id)\
{\
	if (id == -1) {\
		id = pkd->kdTop[iCell].pLower;\
		if (id != -1) iCell = ROOT;\
		else iCell = LOWER(iCell);\
		}\
	else iCell = LOWER(iCell);\
	}

#define pkdUpper(pkd,iCell,id)\
{\
	if (id == -1) {\
		id = pkd->kdTop[iCell].pLower;\
		if (id != -1) iCell = ROOT;\
		else iCell = UPPER(iCell);\
		}\
	else iCell = UPPER(iCell);\
	}

#define pkdParent(pkd,iCell,id)\
{\
	iCell = PARENT(iCell);\
	if (iCell == ROOT) {\
		if (id != -1) {\
			iCell = pkd->piLeaf[id];\
			id = -1;\
			}\
		}\
	}

#define pkdNext(pkd,iCell,id)\
{\
	SETNEXT(iCell);\
	if (iCell == ROOT) {\
		if (id != -1) {\
			iCell = pkd->piLeaf[id];\
			id = -1;\
			SETNEXT(iCell);\
			}\
		}\
	}

double pkdGetTimer(PKD,int);
double pkdGetSystemTimer(PKD,int);
double pkdGetWallClockTimer(PKD,int);
void pkdClearTimer(PKD,int);
void pkdStartTimer(PKD,int);
void pkdStopTimer(PKD,int);
void pkdInitialize(PKD *,MDL,int,int,int,FLOAT *,int,int,int);
void pkdFinish(PKD);
void pkdReadTipsy(PKD,char *,int,int,int,int,double,double);
void pkdSetSoft(PKD pkd,double dSoft);
#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double, double, int);
void pkdPreVariableSoft(PKD pkd,int iVariableSoftType);
void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul,int iVariableSoftType);
#endif
void pkdCalcBound(PKD,BND *,BND *,BND *,BND *);
void pkdGasWeight(PKD);
void pkdRungDDWeight(PKD, int, double);
int pkdWeight(PKD,int,FLOAT,int,int,int,int *,int *,FLOAT *,FLOAT *);
int pkdLowerPart(PKD,int,FLOAT,int,int);
int pkdUpperPart(PKD,int,FLOAT,int,int);
int pkdWeightWrap(PKD,int,FLOAT,FLOAT,int,int,int,int *,int *,FLOAT *,FLOAT *);
int pkdLowerPartWrap(PKD,int,FLOAT,FLOAT,int,int);
int pkdUpperPartWrap(PKD,int,FLOAT,FLOAT,int,int);
int pkdLowerOrdPart(PKD,int,int,int);
int pkdUpperOrdPart(PKD,int,int,int);
int pkdActiveTypeOrder(PKD, unsigned int);
int pkdActiveOrder(PKD);
int pkdColRejects(PKD,int,FLOAT,FLOAT,int);
int pkdSwapRejects(PKD,int);
int pkdSwapSpace(PKD);
int pkdFreeStore(PKD);

// Mike: spends >5% of the time in this due to calls from pst - just make into a macro
//int pkdLocal(PKD);
#define pkdLocal(pkd) ((pkd)->nLocal)

void pkdTotals(PKD pkd, int *nDark, int *nGas, int *nStar);
int pkdActive(PKD);
int pkdTreeActive(PKD);
int pkdInactive(PKD);
int pkdTreeInactive(PKD);
int pkdNodes(PKD);
void pkdDomainColor(PKD);
int pkdColOrdRejects(PKD,int,int);
void pkdLocalOrder(PKD);
void pkdWriteTipsy(PKD,char *,int,int,double,double,int);
void pkdCombine(KDN *,KDN *,KDN *);
void pkdCalcCell(PKD,KDN *,FLOAT *,int,struct pkdCalcCellStruct *);
extern double MPI_Wtime();
double pkdCalcOpen(KDN *,int,double,int);
void pkdBuildLocal(PKD,size_t,int,double,int,int,int,KDN *);
void pkdBuildBinary(PKD,int,int,double,int,int,int,KDN *);
void pkdThreadTree(PKD pkd,int iCell,int iNext);
#ifdef GR_DRAG
void pkdGRDragGetSunAccel(PKD,double,double []);
void pkdGRIntegrateCloseParticles(PKD,double,double,double,int *,double*, double []);
void pkdGRCorrectaSun(PKD, double []);
#endif /*GR_DRAG*/
void pkdGravAll(PKD,int,int,int,int,int,double,double,int,double,double *,int *,
				double *,double *,double *,CASTAT *,double *);
void pkdCalcEandL(PKD,double *,double *,double *,double []);
void pkdCalcEandLExt(PKD,double *,double[],double [],double *);
void pkdDrift(PKD,double,FLOAT *,int,int,FLOAT);
void pkdUpdateUdot(PKD pkd,double,double,double,int,int);
void pkdKick(PKD pkd,double,double, double, double, double, double, int, double, double);
void pkdKickPatch(PKD pkd, double dvFacOne, double dvFacTwo,
#ifdef NEED_VPRED
				  double dvPredFacOne, double dvPredFacTwo,
#endif
				  double dOrbFreq, int bOpen);
void pkdReadCheck(PKD,char *,int,int,int,int);
void pkdWriteCheck(PKD,char *,int,int);
void pkdDistribCells(PKD,int,KDN *);
void pkdCalcRoot(PKD,struct ilCellNewt *);
void pkdDistribRoot(PKD,struct ilCellNewt *);
void pkdSwapAll(PKD pkd, int idSwap);
double pkdMassCheck(PKD pkd);
double pkdGetMassChange(PKD pkd);
void pkdMassMetalsEnergyCheck(PKD pkd, double *dTotMass, double *dTotMetals, 
                    double *dTotOx, double *dTotFe, double *dTotEnergy);
void pkdSetRung(PKD pkd, int iRung);
void pkdBallMax(PKD pkd, int iRung, int bGreater, double ddHonHLimit);
int pkdActiveRung(PKD pkd, int iRung, int bGreater);
int pkdCurrRung(PKD pkd, int iRung);
void pkdGravStep(PKD pkd, double dEta);
void pkdAccelStep(PKD pkd, double dEta, double dVelFac, double
				  dAccFac, int bDoGravity, int bEpsAcc, int bSqrtPhi, double dhMinOverSoft);
void pkdDensityStep(PKD pkd, double dEta, double dRhoFac);
int pkdDtToRung(PKD pkd, int iRung, double dDelta, int iMaxRung, int bAll,
		int *pnMaxRung, int *piMaxRungIdeal);
void pkdInitDt(PKD pkd, double dDelta);
int pkdRungParticles(PKD,int);
void pkdCoolVelocity(PKD,int,double,double,double);
void pkdGrowMass(PKD pkd,int nGrowMass, double dDeltaM);
void pkdInitAccel(PKD);
int pkdOrdWeight(PKD,int,int,int,int,int *,int *);
#ifdef STARFORM
void pkdUnDeleteParticle(PKD pkd, PARTICLE *p);
#endif
void pkdDeleteParticle(PKD pkd, PARTICLE *p);
#ifdef COLLMOD_ZML
void pkdMarkParticleSC(PKD pkd, PARTICLE *p, int an, int bn);
#endif
void pkdNewParticle(PKD pkd, PARTICLE p);
int pkdResetTouchRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdActiveExactType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask);
int pkdActiveType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdSetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdResetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdCountType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask);
int pkdActiveMaskRung(PKD pkd, unsigned int iSetMask, int iRung, int bGreater );
int pkdActiveTypeRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater);
int pkdSetTypeFromFile(PKD pkd, int iSetMask, int biGasOrder, char *file, int *pniOrder, int *pnSet, int *pnSetiGasOrder);

void pkdSetParticleTypes(PKD pkd, int nSuperCool);

struct SoughtParticle {
  int iOrder;
  int iActive; 
  double x,y,z;
};

int pkdSoughtParticleList(PKD pkd, int iTypeSought, int nMax, int *n, struct SoughtParticle *sp);

void pkdCoolUsingParticleList(PKD pkd, int nList, struct SoughtParticle *sp);

void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
				  int *nDeltaStar);
void pkdNewOrder(PKD pkd, int nStart);
void pkdMoveParticle(PKD pkd, double *xcenter,double *xoffset,int iOrder);


struct outGetNParts { 
	int n;
    int nGas;
    int nDark;
    int nStar;
    int iMaxOrderGas;
    int iMaxOrderDark;
    int iMaxOrderStar;
    };

void pkdGetNParts(PKD pkd, struct outGetNParts *out );
void pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar, int, int nMaxOrderGas,
				  int nMaxOrderDark);
void pkdSunIndirect(PKD,double *,int,double,double);
void pkdLogHalo(PKD);
void pkdHernquistSpheroid(PKD pkd);
void pkdNFWSpheroid(PKD pkd, double M_200, double r_200, double c, double dSoft);
void pkdElliptical(PKD pkd, int bEllipticalDarkNFW);
void pkdHomogSpheroid(PKD pkd);
void pkdBodyForce(PKD pkd);
void pkdMiyamotoDisk(PKD pkd);
void pkdTimeVarying(PKD pkd,double dTime);
#ifdef ROT_FRAME
void pkdRotFrame(PKD pkd, double dOmega, double dOmegaDot);
#endif

#ifdef GASOLINE

void pkdUpdateuDot(PKD,double,double,double,int,int);
void pkdUpdateShockTracker(PKD,double, double, double);
void pkdAdiabaticGasPressure(PKD, double gammam1, double gamma);
void pkdCoolingGasPressure(PKD, double gammam1, double gamma);
void pkdLowerSoundSpeed(PKD, double);
void pkdInitEnergy(PKD pkd, double dTuFac, double z, double dTime );
void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta);
int pkdSphCurrRung(PKD pkd, int iRung, int bGreater);
void pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot, int bViscosityLimitdt);
void pkdSphViscosityLimiter(PKD pkd, int bOn, int bShockTracker);

void pkdDensCheck(PKD pkd, int iRung, int bGreater, int iMeasure, void *data);

#endif /* GASOLINE */
#ifdef GLASS
void pkdGlassGasPressure(PKD, void *in);
void pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR);
#endif

#ifdef SLIDING_PATCH
double SHEAR(int,double,PATCH_PARAMS *);
#define SHEAR(ix,t,pp)\
	((ix) < 0 ? fmod(0.5*(pp)->dLength - 1.5*(ix)*(pp)->dOrbFreq*(pp)->dWidth*(t),(pp)->dLength) - 0.5*(pp)->dLength :\
	 (ix) > 0 ? 0.5*(pp)->dLength - fmod(0.5*(pp)->dLength + 1.5*(ix)*(pp)->dOrbFreq*(pp)->dWidth*(t),(pp)->dLength) : 0.0)
void pkdPatch(PKD pkd);
int pkdRandAzWrap(PKD pkd);
#endif

#ifdef COLLISIONS
int pkdNumRejects(PKD pkd);
void pkdReadSS(PKD pkd, char *pszFileName, int nStart, int nLocal);
void pkdWriteSS(PKD pkd, char *pszFileName, int nStart, int bReduced);
void pkdAddUnifGrav(PKD pkd, double dgx, double dgy, double dgz);
void pkdNextEncounter(PKD pkd, double *dt);
void pkdMarkEncounters(PKD pkd, double dt);
#ifdef SIMPLE_GAS_DRAG
void pkdSimpleGasDrag(PKD pkd,int iFlowOpt,int bEpstein,double dGamma,
					  double dTime);
#endif
#ifdef JOHNNY
int pkdFindEscapers(PKD pkd,double dCentMass,double dStarMass);
int pkdCheckForKepler(PKD pkd,double dCentMass,double dDelta);
#endif /* JOHNNY */
#endif /* COLLISIONS */

#ifdef SPRINGS
int pkdBreakSprings(PKD pkd,int iOrder1,int iOrder2);
#endif /* SPRINGS */

#ifdef DEM

void pkdDEMZeroSmallMotions(PKD pkd,double dAccCritSq,double dDeltaSq);
void pkdDEMStats(PKD pkd,DEM_STATS *ds);

#ifdef AGGS
void pkdDEMAggsSetMass(PKD pkd,int iAggIdx,FLOAT fMass);
#endif /* AGGS */

#if defined(WALLS) && defined(WALLS_REACT)
void pkdDEMWallsReact(PKD pkd,double *dTotalZForceFromParticles);
#endif /* defined(WALLS) && defined(WALLS_REACT) */

#endif /* DEM */

#ifdef CHARGE
void pkdChargeInit(PKD pkd);
#endif /* CHARGE */

#ifdef DEM_TIDAL_SPACE
void pkdDEMTidalFindPlanet(PKD pkd,int iColor,int *bFound,double dPos[],double dVel[]);
void pkdDEMTidalFindMarker(PKD pkd,int iOrder,int *bFound,double dPos[],double dVel[]);
#endif /* DEM_TIDAL_SPACE */

#ifdef DEM_TIDAL_LOCAL
void pkdDEMTidal(PKD pkd,const DEM_TIDAL *d0,const DEM_TIDAL *d);
#endif /* DEM_TIDAL_LOCAL */

void pkdMassInR(PKD pkd, double R, double *pdMass, FLOAT *com);

#ifdef NEED_VPRED
#ifdef GASOLINE
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo, double duDelta,
				  int iGasModel, double z, double duDotLimit);
#else
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo);
#endif
#endif

void pkdInitRotBar(PKD pkd, ROTBAR rotbar);
void pkdRotatingBar(PKD pkd, double amp, /* relative amplitude of bar */
		    double posang, /* position angle of bar */
		    double b5,	/* radial scale length (^5) */
		    FLOAT *aCom, /* Center of mass */
		    double *accCom, /* acceleration (returned) */
		    double *dTorque); /* acceleration (returned) */


void pkdCOM(PKD pkd, double *com);
void pkdCOMByType(PKD pkd, int type, double *com);
void pkdOldestStar(PKD pkd, double *com);
int pkdSetSink(PKD pkd, double dSinkMassMin);
#endif

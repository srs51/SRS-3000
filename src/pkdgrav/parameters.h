#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include <sys/param.h> /* for MAXPATHLEN */
#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

#include "cosmo.h"
#include "patch.h"
#include "foverride.h"

#ifdef COLLISIONS
typedef struct {
  int nUnifGrav;
  double *dtUnifGrav,*dxUnifGrav,*dyUnifGrav,*dzUnifGrav;
  } TIME_VAR_PARAMS;
#endif

struct parameters {
	/*
	 ** Parameters for PKDGRAV.
	 */
	int nThreads;
	int bDiag;
	int bOverwrite;
	int bVWarnings;
	int bVStart;
	int bVStep;
	int bVRungStat;
	int bVDetails;
        int bLogTiming;
        int bLogTimingSubStep;
        int bLogTimingStep;
        int bLogTimingSubStepTot;
        int bLogTimingStepTot;
	int bPeriodic;
	int bRestart;
	int bParaRead;
	int bParaWrite;
	int bUseWallClock;
	int bCannonical;
	int bStandard;
	int bKDK;
	int bBinary;
	int bGravStep;
	int bEpsAccStep;
	int bSqrtPhiStep;
	int bAccelStep; /* true if bEpsAccStep or bSqrtPhiStep */
	int bDensityStep;
    int bDeltaAccelStep;
    int bDeltaAccelStepGasTree;
	int nTruncateRung;
	int bNonSymp;
    int iBinaryOutput;
    int bPackedVector;
	int bDoDensity;
 	int iReadIOrder;
 	int bDoIOrderOutput;
 	int bDohOutput;
 	int bDoSphhOutput;
	int bDodtOutput;
	int bDoIonOutput;
	int bDoDminOutput;
	int bSymCool;
        int bDoGravity;
        int bDoSelfGravity;
        int bFandG;
        int bHeliocentric;
            double dSunSoft;
        int bLogHalo;
        int bHernquistSpheroid;
        int bNFWSpheroid;
            double dNFWm200;
            double dNFWr200;
            double dNFWconc;
            double dNFWsoft;
        int bElliptical;
        int bEllipticalDarkNFW;
            int bHomogSpheroid;
        int bBodyForce;
	int bMiyamotoDisk;
	int bTimeVarying;
	int bRotatingBar;
        ROTBAR  rotbar;
#ifdef ROT_FRAME
	int bRotFrame;
	double dOmega;
	double dOmegaDot;
#endif
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
#ifdef SIMPLE_GAS_DRAG
	int bSimpleGasDrag;
	int bEpstein;
	double dGamma;
#endif

	int nBucket;
	int iOutInterval;
#ifdef COLLISIONS
  int iRedOutInterval;
#endif
	int iLogInterval;
	int iCheckInterval;
	int iOrder;
	int bEwald;
	int iEwOrder;
	int nReplicas;
	int iStartStep;
	int nSteps;
	int bCollDelay;
	double dCollTimer;
	int nSmooth;
	int iMaxRung;
	int nSuperCool;
	int nGrowMass;
	int iWallRunTime;
	int bPhysicalSoft;  
	int bSoftMaxMul;
	int bVariableSoft;
	int nSoftNbr;
	int bSoftByType;
    int bVariableSoftStar;
    int bVariableSoftGas;
    int bVariableSoftDark;
	int bDoSoftOutput;
    int bDoSinks;
    int bBHSink;
    int bDoSinksAtStart;
    int bSinkThermal;
    int iSinkRung;
	double dEta;
        double dEtaDeltaAccel;
	double dExtraStore;
	double dSoft;
	double dSoftMax;
	double dDelta;
	double dEwCut;
	double dEwhCut;
	double dTheta;
	double dTheta2;
	double daSwitchTheta;
	double dAbsPartial;
	double dRelPartial;
	double dAbsTotal;
	double dRelTotal;
	double dPeriod;
	double dxPeriod;
	double dyPeriod;
	double dzPeriod;
	CSM csm;
	double dRedTo;
	double dCentMass;
    double dSinkRadius;
    double dSinkBoundOrbitRadius;
    double dDeltaSink;
    double dSinkMassMin;
    double dBHSinkEddEff;
    double dBHSinkFeedbackEff;
    double dBHSinkEddFactor;
    double dBHSinkFeedbackFactor;
    double dBHSinkAlpha;
    double dSinkCurrentDelta;
	char achDigitMask[MAXPATHLEN];
	char achInFile[MAXPATHLEN];
	char achOutName[MAXPATHLEN];
	char achDataSubPath[MAXPATHLEN];
	double dCoolFac;
	double dCoolDens;
	double dCoolMaxDens;
	double dGrowDeltaM;
	double dGrowStartT;
	double dGrowEndT;
	double dFracNoDomainDecomp;
	double dFracNoDomainDimChoice;
	int    bRungDD;
	double dRungDDWeight;
	/*
	 ** Additional parameters for GASOLINE.
	 */
	int bGeometric;
	int bGasAdiabatic;
	int bGasIsothermal;
	int bGasCooling;
	int bGasCoolingNonEqm;
	int iGasModel;
	double dEtaCourant;
	double dEtauDot;
	double duDotLimit;
	double dShockTrackerA;
	double dShockTrackerB;
	double dConstAlpha;
	double dConstBeta;
	double dConstGamma;
	double dMeanMolWeight;
	double dGasConst;
	double dMsolUnit;
	double dKpcUnit;
	double ddHonHLimit;
	double dGmPerCcUnit;
	double dComovingGmPerCcUnit;
	double dErgPerGmUnit;
	double dSecUnit;
	int    bViscosityLimiter;
	int    bViscosityLimitdt;
	int    bShockTracker;
	int    bBulkViscosity;
	int    bGasDomainDecomp;
	int    bLowerSoundSpeed;
	int    bFastGas;
	double dFracFastGas;
	double dhMinOverSoft;
	int    bDoGas;
	int    bSphStep;
#if defined(GASOLINE) && !defined(NOCOOLING)
	COOLPARAM CoolParam;
#endif
	int    bSN;
	double dSNRhoCut;
 	double dSNTMin;
	double dSNTMax;
	double dSNMetalCut;
	double dSNHeatFraction;
	double dDumpFrameStep;
	double dDumpFrameTime;
	int    bStarForm;
	int    bFeedBack;
	int    bFormOutputs;
#ifdef SIMPLESF
	double SSF_dEfficiency;
    double SSF_dTMax;
    double SSF_dPhysDenMin;
    double SSF_dComovingDenMin;
    double SSF_dESNPerStarMass;
    double SSF_dInitStarMass;
	double SSF_dtCoolingShutoff;
    int SSF_bdivv;
#endif
#ifdef STARFORM
	STFM   stfm;
	FB     fb;
        SN sn;
        double dDeltaStarForm;
        int bSNTurnOffCooling;
	int bShortCoolShutoff;
	int bSmallSNSmooth;
        int iStarFormRung;
#endif
	double dKBoltzUnit;
#ifdef GLASS
	/*
	 ** Additional parameters for GLASS.
	 */
	double dGlassDamper;
	/* Hack for shock tube glasses */
	double dGlassPoverRhoL;
	double dGlassPoverRhoR;
	double dGlassxL;
	double dGlassxR;
	double dGlassVL;
	double dGlassVR;
#endif
#ifdef COLLISIONS
  /*
  ** Additional parameters for collision code...
  */
  int bAllowSimulColl;
  int bFindRejects;
  int iCollLogOption;
  int iMinCollRung;
  char achCollLog[MAXPATHLEN];
  double dSmallStep;
  double dxUnifGrav;
  double dyUnifGrav;
  double dzUnifGrav;
  int    iMinBinaryRung;
  double dBallVelFact;
  double dMaxBinaryEcc;
  TIME_VAR_PARAMS sTimeVarParams;
#ifdef SLIDING_PATCH
  int    iRandStep;
  double dLargeMass;
  double dRandBall;
  int iNextRandomization;
#endif
  COLLISION_PARAMS CP;
#endif /* COLLISIONS */
#ifdef SPECIAL_PARTICLES
  int nSpecial;
  int iSpecialId[MAX_NUM_SPECIAL_PARTICLES];
  SPECIAL_PARTICLE_DATA sSpecialData[MAX_NUM_SPECIAL_PARTICLES];
#endif /* SPECIAL_PARTICLES */
  /*
  ** Additional parameters for force override options...
  */
  FOVERRIDE_PARAMS FO;
  };

#endif

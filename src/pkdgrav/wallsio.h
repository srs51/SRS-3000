#ifndef WALLSIO_HINCLUDED
#define WALLSIO_HINCLUDED

#define WALLS_MAX_STR_LEN 4096 /* must be at least 1 */

/* wall types */

enum {WallPlane=0,WallTriangle,WallRectangle,WallDisk,WallCylinderInfinite,WallCylinderFinite,WallShell};

/* verbosity flags */

#define WALLS_SILENT  0
#define WALLS_VERBOSE 1

/* wall data struct */

typedef struct {
  int iType;
  double vOrigin[3]; /* can treat these as 3-vectors */
  double vOrient[3];
  double vVertex1[3];
  double vVertex2[3];
  double vVel[3];
  double dOscAmp;
  double dOscFreq;
  double vOscVec[3];
  double dRadius;
  double dHoleRadius;
  double dLength;
  double dTaper;
  double dOpenAngle;
  double dAngSpeed;
  double dEpsN;
  double dEpsT;
  double dKn;
  double dKt;
  double dMuS;
  double dMuR;
  double dMuT;
  double dKnOuter;
  double dKtOuter;
  double dKnInner;
  double dKtInner;
  double dInnerOverlapBoundary;
  int iColor;
  double dTrans;
  double dMass;
  } WALL_DATA;

/* function prototype */

int wallsParseWallsFile(FILE *fp,int *nWalls,WALL_DATA **pWallsData,double *dTime,int bVerbose);

#endif /* !WALLSIO_HINCLUDED */

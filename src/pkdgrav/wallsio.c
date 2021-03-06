/*
** wallsio.c -- DCR 7/14/10
** =========
** Parses wall data for simulation (pkdgrav) and visualization (ssdraw).
** Documentation available in walls.pdf (standard ss_core distribution).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <assert.h>
#include "wallsio.h"

/*#define STANDALONE*/ /* uncomment for standalone testing */

/* (to compile standalone: gcc -I$ss_dir/include $ss_dir/wallsio.c -lm) */

#define TOLERANCE 0.001 /* max rectangle vertex vectors dot product mag */

/* initializations for tokens */

enum {TokenNot=-1,TokenTime,TokenLengthUnit,TokenMassUnit,TokenTimeUnit,TokenDefaults,TokenWall,TokenType,TokenOrigin,TokenOrient,TokenVertex1,TokenVertex2,TokenVel,TokenOscAmp,TokenOscFreq,TokenOscVec,TokenRadius,TokenHoleRadius,TokenLength,TokenTaper,TokenOpenAngle,TokenAngSpeed,TokenEpsN,TokenEpsT,TokenKn,TokenKt,TokenKnOuter,TokenKtOuter,TokenKnInner,TokenKtInner,TokenInnerOverlapBoundary,TokenMuS,TokenMuR,TokenMuT,TokenColor,TokenTrans,TokenMass}; /* must match entries in achToken (except TokenNot) */

#define NUM_TOKENS 36 /* must match number of entries in achToken (except TokenNot) */

static const char *achToken[] = {"time","lengthunit","massunit","timeunit","defaults","wall","type","origin","orient","vertex1","vertex2","velocity","osc-ampl","osc-freq","osc-vec","radius","hole-radius","length","taper","open-angle","ang-speed","epsn","epst","k_n","k_t","k_n_outer","k_t_outer","k_n_inner","k_t_inner","inner_overlap_boundary","mu_s","mu_r","mu_t","color","transparency","mass"};

/* linear algebra stuff */

typedef double Vec[3];

static void vecSet(Vec v,double x,double y,double z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
	}

static void vecZero(Vec v)
{
	vecSet(v,0.0,0.0,0.0);
	}

static void vecCopy(const Vec vSrc,Vec vDst)
{
	vDst[0] = vSrc[0];
	vDst[1] = vSrc[1];
	vDst[2] = vSrc[2];
	}

static void vecScale(Vec v,double dScale)
{
	v[0] *= dScale;
	v[1] *= dScale;
	v[2] *= dScale;
	}

/*DEBUG not used...
static void vecSub(const Vec v1,const Vec v2,Vec v)
{
	v[0] = v1[0] - v2[0];
	v[1] = v1[1] - v2[1];
	v[2] = v1[2] - v2[2];
	}
*/

static double vecDot(const Vec v1,const Vec v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
	}

static double vecMagSq(const Vec v)
{
	return vecDot(v,v);
	}

static double vecMag(const Vec v)
{
	return sqrt(vecMagSq(v));
	}

static void vecNorm(Vec v)
{
	double dMag = vecMag(v);

	assert(dMag > 0.0);
	vecScale(v,1.0/dMag);
	}

static void vecCross(const Vec v1,const Vec v2,Vec v)
{
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

/* initializations for types */

enum {TypeNot=-1,TypePlane=WallPlane,TypeTriangle,TypeRectangle,TypeDisk,TypeCylinderInfinite,TypeCylinderFinite,TypeShell}; /* must match entries in achType (except TypeNot) */

#define NUM_TYPES 7 /* must match number of entries in achType */

static const char *achType[] = {"plane","triangle","rectangle","disk","cylinder-infinite","cylinder-finite","shell"};

static char *nextNonWhitespace(char *p)
{
	while (isspace(*p) || *p == ',' || *p == '=') /* matches space, tab, newline, comma, equals */
		++p;

	return p;
	}

static char *nextWhitespace(char *p)
{
	while (*p != '\0' && (!isspace(*p) && *p != ',' && *p != '='))
		++p;

	return p;
	}

static char *getNextWord(char *p,char **pWord)
{
	*pWord = nextNonWhitespace(p);

	if (**pWord == '\0' || **pWord == '#' || **pWord == '!') /* EOL or comment */
		return NULL; /* no word found */

	p = nextWhitespace(*pWord);

	*p = '\0';

	return p + 1; /* ok because '\n' guaranteed to be at EOL */
	}

static char *getVec(char *p,Vec vVec)
{
	char *pWord;
	int k;

	for (k=0;k<3;k++) {
		if ((p = getNextWord(p,&pWord)) == NULL)
			return NULL;
		vVec[k] = atof(pWord);
		}

	return p;
	}

static int getToken(const char *achWord)
{
	int i;

	for (i=0;i<NUM_TOKENS;i++)
		if (strcasecmp(achWord,achToken[i]) == 0)
			return i;

	return TokenNot; /* no match */
	}

static int getType(const char *achWord)
{
	int i;

	for (i=0;i<NUM_TYPES;i++)
		if (strcasecmp(achWord,achType[i]) == 0)
			return i;

	return TypeNot; /* no match */
	}

static void parseErrorEOL(int iLine,int iToken)
{
	fprintf(stderr,"Line %i: missing argument for token \"%s\" (end of line encountered instead).\n",iLine,achToken[iToken]);
	}

static void parseErrorArg(int iLine,const char *achArg,int iToken)
{
	fprintf(stderr,"Line %i: invalid argument ",iLine);
	if (achArg != NULL)
		fprintf(stderr,"(\"%s\") ",achArg);
	fprintf(stderr,"for token \"%s\".\n",achToken[iToken]);
	}

static void parseErrorMode(int iLine,int iToken)
{
	fprintf(stderr,"Line %i: token \"%s\" not in wall or defaults block.\n",iLine,achToken[iToken]);
	}

static void parseWarnType(int iLine,int iType,int iToken)
{
	fprintf(stderr,"Line %i: WARNING: token \"%s\" not used for this wall type (\"%s\").\n",iLine,achToken[iToken],achType[iType]);
	}

int wallsParseWallsFile(FILE *fp,int *nWalls,WALL_DATA **pWallData,double *dTime,int bVerbose)
{
	enum {ModeNone=-1,ModeDefaults,ModeWall};

	WALL_DATA *pWall=NULL,wDflt;
	Vec vTmp;
	double dLengthUnit,dMassUnit,dTimeUnit,dTmp;
	char achLine[WALLS_MAX_STR_LEN],*pWord,*p;
	int i,iLine=0,iMode=ModeNone,iTmp;

	*nWalls = 0;
	*pWallData = NULL;
	*dTime = 0.0;

	dLengthUnit = dMassUnit = dTimeUnit = 1.0;

	wDflt.iType = TypePlane;
	vecZero(wDflt.vOrigin);
	vecSet(wDflt.vOrient,0.0,0.0,1.0); /* must be unit vector */
	vecSet(wDflt.vVertex1,1.0,0.0,0.0);
	vecSet(wDflt.vVertex2,0.0,1.0,0.0);
	vecZero(wDflt.vVel);
	wDflt.dOscAmp = 0.0;
	wDflt.dOscFreq = 0.0;
	vecCopy(wDflt.vOrient,wDflt.vOscVec); /* must be unit vector */
	wDflt.dRadius = 1.0;
	wDflt.dHoleRadius = 0.0;
	wDflt.dLength = 1.0;
	wDflt.dTaper = 0.0;
	wDflt.dOpenAngle = 0.0;
	wDflt.dAngSpeed = 0.0;
	wDflt.dEpsN = 1.0; /* note these override particle values! */
	wDflt.dEpsT = 1.0;
	/*
	** The next ten parameters are for DEM.  The latter five are for when the DEM_TWOLAYERS macro is defined.
    ** The defaults indicate that the corresponding parameters inherit particle values.
	*/
	wDflt.dKn = 0.0;
	wDflt.dKt = 0.0;
	wDflt.dMuS = -1.0;
	wDflt.dMuR = -1.0;
	wDflt.dMuT = -1.0;
	wDflt.dKnOuter = 0.0;
	wDflt.dKtOuter = 0.0;
	wDflt.dKnInner = 0.0;
	wDflt.dKtInner = 0.0;
	wDflt.dInnerOverlapBoundary = -1.0;
	/* remaining parameters */
	wDflt.iColor = 223; /* light gray */
	wDflt.dTrans = 0.0; /* opaque */
	wDflt.dMass = 0.0;

	while (fgets(achLine,WALLS_MAX_STR_LEN,fp) != NULL) {

		++iLine;

		/* some sanity checks */

		i = strlen(achLine);
		if (i == 0)
			continue; /* nothing on this line -- skip it! */
		if (i == WALLS_MAX_STR_LEN && achLine[i - 1] != '\n') {
			fprintf(stderr,"Line %i: maximum number of allowed characters exceeded (%i).\n",iLine,WALLS_MAX_STR_LEN - 1); /* -1 to allow for carriage return */
			return 1;
			}
		assert(achLine[i - 1] == '\n'); /* otherwise fgets() not working! */

		p = achLine;
		while ((p = getNextWord(p,&pWord)) != NULL) {
			switch (getToken(pWord)) {
			case TokenTime:
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenTime);
					return 1;
					}
				*dTime = atof(pWord)*dTimeUnit;
				if (iMode != ModeDefaults) /* allow this token to be set in defaults mode */
					iMode = ModeNone;
				pWall = NULL;
				break;
			case TokenLengthUnit:
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenLengthUnit);
					return 1;
					}
				dLengthUnit = atof(pWord);
				if (dLengthUnit <= 0.0) {
					parseErrorArg(iLine,pWord,TokenLengthUnit);
					return 1;
					}
				if (iMode != ModeDefaults)
					iMode = ModeNone;
				pWall = NULL;
				break;
			case TokenMassUnit:
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenMassUnit);
					return 1;
					}
				dMassUnit = atof(pWord);
				if (dMassUnit <= 0.0) {
					parseErrorArg(iLine,pWord,TokenMassUnit);
					return 1;
					}
				if (iMode != ModeDefaults)
					iMode = ModeNone;
				pWall = NULL;
				break;
			case TokenTimeUnit:
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenTimeUnit);
					return 1;
					}
				dTimeUnit = atof(pWord);
				if (dTimeUnit <= 0.0) {
					parseErrorArg(iLine,pWord,TokenTimeUnit);
					return 1;
					}
				if (iMode != ModeDefaults)
					iMode = ModeNone;
				pWall = NULL;
				break;
			case TokenDefaults:
				iMode = ModeDefaults;
				pWall = NULL;
				break;
			case TokenWall:
				++(*nWalls);
				*pWallData = realloc(*pWallData,(*nWalls)*sizeof(WALL_DATA));
				assert(*pWallData != NULL);
				pWall = &((*pWallData)[*nWalls - 1]);
				*pWall = wDflt; /* struct copy */
				iMode = ModeWall;
				break;
			case TokenType:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenType);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenType);
					return 1;
					}
				iTmp = getType(pWord);
				if (iTmp == TypeNot) {
					fprintf(stderr,"Line %i: unrecognized wall type (\"%s\").\n",iLine,pWord);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.iType = iTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->iType = iTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenOrigin:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenOrigin);
					return 1;
					}
				if ((p = getVec(p,vTmp)) == NULL) {
					parseErrorEOL(iLine,TokenOrigin);
					return 1;
					}
				vecScale(vTmp,dLengthUnit);
				switch (iMode) {
				case ModeDefaults:
					vecCopy(vTmp,wDflt.vOrigin);
					break;
				case ModeWall:
					assert(pWall != NULL);
					vecCopy(vTmp,pWall->vOrigin);
					break;
				default:
					assert(0);
					}
				break;
			case TokenOrient:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenOrient);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypePlane && pWall->iType != TypeDisk && pWall->iType != TypeCylinderInfinite && pWall->iType != TypeCylinderFinite && pWall->iType != TypeShell)
						parseWarnType(iLine,pWall->iType,TokenOrient);
					}
				if ((p = getVec(p,vTmp)) == NULL) {
					parseErrorEOL(iLine,TokenOrient);
					return 1;
					}
				dTmp = vecMag(vTmp);
				if (dTmp == 0.0) {
					parseErrorArg(iLine,NULL,TokenOrient);
					return 1;
					}
				if (dTmp != 1.0) {
					vecNorm(vTmp);
					if (bVerbose)
						fprintf(stderr,"Line %i: WARNING: normalized %s to unit vector.\n",iLine,achToken[TokenOrient]);
					}
				switch (iMode) {
				case ModeDefaults:
					vecCopy(vTmp,wDflt.vOrient);
					break;
				case ModeWall:
					assert(pWall != NULL);
					vecCopy(vTmp,pWall->vOrient);
					break;
				default:
					assert(0);
					}
				break;
			case TokenVertex1:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenVertex1);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeTriangle && pWall->iType != TypeRectangle)
						parseWarnType(iLine,pWall->iType,TokenVertex1);
					}
				if ((p = getVec(p,vTmp)) == NULL) {
					parseErrorEOL(iLine,TokenVertex1);
					return 1;
					}
				if (vecMag(vTmp) == 0.0) {
					parseErrorArg(iLine,NULL,TokenVertex1);
					return 1;
					}
				vecScale(vTmp,dLengthUnit);
				switch (iMode) {
				case ModeDefaults:
					vecCopy(vTmp,wDflt.vVertex1);
					break;
				case ModeWall:
					assert(pWall != NULL);
					vecCopy(vTmp,pWall->vVertex1);
					break;
				default:
					assert(0);
					}
				break;
			case TokenVertex2:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenVertex2);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeTriangle && pWall->iType != TypeRectangle)
						parseWarnType(iLine,pWall->iType,TokenVertex2);
					}
				if ((p = getVec(p,vTmp)) == NULL) {
					parseErrorEOL(iLine,TokenVertex2);
					return 1;
					}
				if (vecMag(vTmp) == 0.0) {
					parseErrorArg(iLine,NULL,TokenVertex2);
					return 1;
					}
				vecScale(vTmp,dLengthUnit);
				switch (iMode) {
				case ModeDefaults:
					vecCopy(vTmp,wDflt.vVertex2);
					break;
				case ModeWall:
					assert(pWall != NULL);
					vecCopy(vTmp,pWall->vVertex2);
					break;
				default:
					assert(0);
					}
				break;
			case TokenVel:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenVel);
					return 1;
					}
				if ((p = getVec(p,vTmp)) == NULL) {
					parseErrorEOL(iLine,TokenVel);
					return 1;
					}
				vecScale(vTmp,dLengthUnit/dTimeUnit);
				switch (iMode) {
				case ModeDefaults:
					vecCopy(vTmp,wDflt.vVel);
					break;
				case ModeWall:
					assert(pWall != NULL);
					vecCopy(vTmp,pWall->vVel);
					break;
				default:
					assert(0);
					}
				break;
			case TokenOscAmp:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenOscAmp);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenOscAmp);
					return 1;
					}
				dTmp = atof(pWord)*dLengthUnit;
				switch (iMode) {
				case ModeDefaults:
					wDflt.dOscAmp = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dOscAmp = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenOscFreq:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenOscFreq);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenOscFreq);
					return 1;
					}
				dTmp = atof(pWord)/dTimeUnit;
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenOscFreq);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dOscFreq = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dOscFreq = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenOscVec:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenOscVec);
					return 1;
					}
				if ((p = getVec(p,vTmp)) == NULL) {
					parseErrorEOL(iLine,TokenOscVec);
					return 1;
					}
				dTmp = vecMag(vTmp);
				if (dTmp == 0.0) {
					parseErrorArg(iLine,NULL,TokenOscVec);
					return 1;
					}
				if (dTmp != 1.0) {
					vecNorm(vTmp);
					if (bVerbose)
						fprintf(stderr,"Line %i: WARNING: normalized %s to unit vector.\n",iLine,achToken[TokenOscVec]);
					}
				switch (iMode) {
				case ModeDefaults:
					vecCopy(vTmp,wDflt.vOscVec);
					break;
				case ModeWall:
					assert(pWall != NULL);
					vecCopy(vTmp,pWall->vOscVec);
					break;
				default:
					assert(0);
					}
				break;
			case TokenRadius:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenRadius);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeDisk && pWall->iType != TypeCylinderInfinite && pWall->iType != TypeCylinderFinite && pWall->iType != TypeShell)
						parseWarnType(iLine,pWall->iType,TokenRadius);
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenRadius);
					return 1;
					}
				dTmp = atof(pWord)*dLengthUnit;
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenRadius);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dRadius = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dRadius = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenHoleRadius:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenHoleRadius);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeDisk)
						parseWarnType(iLine,pWall->iType,TokenHoleRadius);
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenHoleRadius);
					return 1;
					}
				dTmp = atof(pWord)*dLengthUnit;
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenHoleRadius);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dHoleRadius = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dHoleRadius = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenLength:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenLength);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeCylinderFinite)
						parseWarnType(iLine,pWall->iType,TokenLength);
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenLength);
					return 1;
					}
				dTmp = atof(pWord)*dLengthUnit;
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenLength);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dLength = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dLength = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenTaper:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenTaper);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeCylinderFinite)
						parseWarnType(iLine,pWall->iType,TokenTaper);
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenTaper);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 || dTmp > 1.0) {
					parseErrorArg(iLine,pWord,TokenTaper);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dTaper = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dTaper = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenOpenAngle:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenOpenAngle);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypeShell)
						parseWarnType(iLine,pWall->iType,TokenOpenAngle);
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenOpenAngle);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 || dTmp > 180.0) {
					parseErrorArg(iLine,pWord,TokenOpenAngle);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dOpenAngle = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dOpenAngle = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenAngSpeed:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenAngSpeed);
					return 1;
					}
				if (iMode == ModeWall) {
					assert(pWall != NULL);
					if (bVerbose && pWall->iType != TypePlane && pWall->iType != TypeDisk && pWall->iType != TypeCylinderInfinite && pWall->iType != TypeCylinderFinite && pWall->iType != TypeShell)
						parseWarnType(iLine,pWall->iType,TokenAngSpeed);
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenAngSpeed);
					return 1;
					}
				dTmp = atof(pWord)/dTimeUnit;
				switch (iMode) {
				case ModeDefaults:
					wDflt.dAngSpeed = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dAngSpeed = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenEpsN:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenEpsN);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenEpsN);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp > 1.0) {
					parseErrorArg(iLine,pWord,TokenEpsN);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dEpsN = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dEpsN = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenEpsT:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenEpsT);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenEpsT);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < -1.0 || dTmp > 1.0) {
					parseErrorArg(iLine,pWord,TokenEpsT);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dEpsT = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dEpsT = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenKn:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenKn);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenKn);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenKn);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dKn = dTmp;
					wDflt.dKt = (2.0/7.0)*wDflt.dKn; /* k_t defaults to 2/7 k_n */
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dKn = dTmp;
					pWall->dKt = (2.0/7.0)*pWall->dKn; /* k_t defaults to 2/7 k_n */
					break;
				default:
					assert(0);
					}
				break;
			case TokenKt:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenKt);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenKt);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenKt);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dKt = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dKt = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenKnOuter:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenKnOuter);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenKnOuter);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenKnOuter);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dKnOuter = dTmp;
					wDflt.dKtOuter = (2.0/7.0)*wDflt.dKnOuter; /* k_t_outer defaults to 2/7 k_n_outer */
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dKnOuter = dTmp;
					pWall->dKtOuter = (2.0/7.0)*pWall->dKnOuter; /* k_t_outer defaults to 2/7 k_n_outer */
					break;
				default:
					assert(0);
					}
				break;
			case TokenKtOuter:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenKtOuter);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenKtOuter);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenKtOuter);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dKtOuter = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dKtOuter = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenKnInner:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenKnInner);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenKnInner);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenKnInner);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dKnInner = dTmp;
					wDflt.dKtInner = (2.0/7.0)*wDflt.dKnInner; /* k_t_inner defaults to 2/7 k_n_inner */
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dKnInner = dTmp;
					pWall->dKtInner = (2.0/7.0)*pWall->dKnInner; /* k_t_inner defaults to 2/7 k_n_inner */
					break;
				default:
					assert(0);
					}
				break;
			case TokenKtInner:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenKtInner);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenKtInner);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenKtInner);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dKtInner = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dKtInner = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenInnerOverlapBoundary:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenInnerOverlapBoundary);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenInnerOverlapBoundary);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 && dTmp != -1.0) {
					parseErrorArg(iLine,pWord,TokenInnerOverlapBoundary);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dInnerOverlapBoundary = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dInnerOverlapBoundary = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenMuS:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenMuS);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenMuS);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 && dTmp != -1.0) {
					parseErrorArg(iLine,pWord,TokenMuS);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dMuS = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dMuS = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenMuR:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenMuR);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenMuR);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 && dTmp != -1.0) {
					parseErrorArg(iLine,pWord,TokenMuR);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dMuR = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dMuR = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenMuT:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenMuT);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenMuT);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 && dTmp != -1.0) {
					parseErrorArg(iLine,pWord,TokenMuT);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dMuT = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dMuT = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenColor:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenColor);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenColor);
					return 1;
					}
				iTmp = atoi(pWord);
				if (iTmp < 0 || iTmp > 255) {
					parseErrorArg(iLine,pWord,TokenColor);
					return 1;
					}
				if ((iTmp == 0 || iTmp == 16) && bVerbose)
					fprintf(stderr,"Line %i: WARNING: drawing color black.\n",iLine);
				switch (iMode) {
				case ModeDefaults:
					wDflt.iColor = iTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->iColor = iTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenTrans:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenTrans);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenTrans);
					return 1;
					}
				dTmp = atof(pWord);
				if (dTmp < 0.0 || dTmp > 100.0) {
					parseErrorArg(iLine,pWord,TokenTrans);
					return 1;
					}
				dTmp = (dTmp > 1.0 ? 0.01*dTmp : dTmp); /* rescale */
				switch (iMode) {
				case ModeDefaults:
					wDflt.dTrans = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dTrans = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenMass:
				if (iMode != ModeDefaults && iMode != ModeWall) {
					parseErrorMode(iLine,TokenMass);
					return 1;
					}
				if ((p = getNextWord(p,&pWord)) == NULL) {
					parseErrorEOL(iLine,TokenMass);
					return 1;
					}
				dTmp = atof(pWord)*dMassUnit;
				if (dTmp < 0.0) {
					parseErrorArg(iLine,pWord,TokenMass);
					return 1;
					}
				switch (iMode) {
				case ModeDefaults:
					wDflt.dMass = dTmp;
					break;
				case ModeWall:
					assert(pWall != NULL);
					pWall->dMass = dTmp;
					break;
				default:
					assert(0);
					}
				break;
			case TokenNot:
				fprintf(stderr,"Line %i: Missing token (found \"%s\" instead).\n",iLine,pWord);
				return 1;
			default:
				assert(0);
				}
			} /* while (words) */
		} /* while (lines) */

	/* final error check, etc. */

	for (i=0;i<*nWalls;i++) {
		pWall = &((*pWallData)[i]);
		if (pWall->iType == TypeTriangle || pWall->iType == TypeRectangle) { /*DEBUG! allowing triangles for soft sphere experimentation*/
			if (pWall->iType == TypeRectangle)
				if (fabs(vecDot(pWall->vVertex1,pWall->vVertex2)) > TOLERANCE) {
					fprintf(stderr,"Wall %i: rectangle axis vectors not perpendicular (|%g| should be <= %g).\n",i,vecDot(pWall->vVertex1,pWall->vVertex2),TOLERANCE);
					return 1;
				}
			/* construct orientation (unit) vector for rectangle */
			vecCross(pWall->vVertex1,pWall->vVertex2,pWall->vOrient);
			dTmp = vecMag(pWall->vOrient);
			assert(dTmp > 0.0); /* error condition already checked for */
			vecScale(pWall->vOrient,1.0/dTmp); /* equivalent to vecNorm(pWall->vOrient) */
		}
		if (pWall->iType == TypeDisk && ((pWall->dRadius > 0.0 && pWall->dHoleRadius >= pWall->dRadius) || (pWall->dRadius == 0.0 && pWall->dHoleRadius > 0.0))) {
			fprintf(stderr,"Wall %i: disk radius must be larger than hole radius.\n",i);
			return 1;
			}
		if (pWall->iType == TypeCylinderFinite && pWall->dLength == 0.0 && pWall->dTaper > 0.0) {
			fprintf(stderr,"Wall %i: zero-length cylinders must have zero taper.\n",i); /* maybe could allow it, but why bother? */
			return 1;
			}
		if (bVerbose && pWall->dOscAmp != 0.0 && pWall->dOscFreq == 0.0)
			fprintf(stderr,"Wall %i: WARNING: zero oscillation frequency (amplitude %g).\n",i,pWall->dOscAmp);
	}
	return 0;
}

#ifdef STANDALONE
int main(int argc,char *argv[1])
{
	WALL_DATA *pWallData;
	double dTime;
	int nWalls;

	FILE *fp;

	assert(WALLS_MAX_STR_LEN > 0);

	if (argc != 2) {
		fprintf(stderr,"Usage: %s walls-file\n",argv[0]);
		return 1;
		}

	printf("Reading wall data from \"%s\"...\n",argv[1]);

	if ((fp = fopen(argv[1],"r")) == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for reading.\n",argv[1]);
		return 1;
		}

	if (wallsParseWallsFile(fp,&nWalls,&pWallData,&dTime,WALLS_VERBOSE) != 0) {
		fprintf(stderr,"Error occurred while parsing walls data.\n");
		return 1;
		}

	fclose(fp);

	printf("Number of walls read = %i\n",nWalls);

	if (nWalls > 0) {
		WALL_DATA *w;
		int i;
		printf("Start time = %g\n",dTime);
		for (i=0;i<nWalls;i++) {
			w = &pWallData[i];
			printf("*** WALL %i ***\n",i);
			printf("%s = %s\n",achToken[TokenType],achType[w->iType]);
			printf("%s = %g,%g,%g\n",achToken[TokenOrigin],w->vOrigin[0],w->vOrigin[1],w->vOrigin[2]);
			printf("%s = %g,%g,%g",achToken[TokenOrient],w->vOrient[0],w->vOrient[1],w->vOrient[2]);
			if (w->iType == WallTriangle || w->iType == WallRectangle)
				printf(" (by construction)");
			printf("\n");
			printf("%s = %g,%g,%g\n",achToken[TokenVertex1],w->vVertex1[0],w->vVertex1[1],w->vVertex1[2]);
			printf("%s = %g,%g,%g\n",achToken[TokenVertex2],w->vVertex2[0],w->vVertex2[1],w->vVertex2[2]);
			printf("%s = %g,%g,%g\n",achToken[TokenVel],w->vVel[0],w->vVel[1],w->vVel[2]);
			printf("%s = %g\n",achToken[TokenOscAmp],w->dOscAmp);
			printf("%s = %g\n",achToken[TokenOscFreq],w->dOscFreq);
			printf("%s = %g,%g,%g\n",achToken[TokenOscVec],w->vOscVec[0],w->vOscVec[1],w->vOscVec[2]);
			printf("%s = %g\n",achToken[TokenRadius],w->dRadius);
			printf("%s = %g\n",achToken[TokenHoleRadius],w->dHoleRadius);
			printf("%s = %g\n",achToken[TokenLength],w->dLength);
			printf("%s = %g\n",achToken[TokenTaper],w->dTaper);
			printf("%s = %g\n",achToken[TokenOpenAngle],w->dOpenAngle);
			printf("%s = %g\n",achToken[TokenAngSpeed],w->dAngSpeed);
			printf("%s = %g\n",achToken[TokenEpsN],w->dEpsN);
			printf("%s = %g\n",achToken[TokenEpsT],w->dEpsT);
			printf("%s = %g\n",achToken[TokenKn],w->dKn);
			printf("%s = %g\n",achToken[TokenKt],w->dKt);
			printf("%s = %g\n",achToken[TokenKnOuter],w->dKnOuter);
			printf("%s = %g\n",achToken[TokenKtOuter],w->dKtOuter);
			printf("%s = %g\n",achToken[TokenKnInner],w->dKnInner);
			printf("%s = %g\n",achToken[TokenKtInner],w->dKtInner);
			printf("%s = %g\n",achToken[TokenInnerOverlapBoundary],w->dInnerOverlapBoundary);
			printf("%s = %g\n",achToken[TokenMuS],w->dMuS);
			printf("%s = %g\n",achToken[TokenMuR],w->dMuR);
			printf("%s = %g\n",achToken[TokenMuT],w->dMuT);
			printf("%s = %i\n",achToken[TokenColor],w->iColor);
			printf("%s = %g\n",achToken[TokenTrans],w->dTrans);
			printf("%s = %g\n",achToken[TokenMass],w->dMass);
			}
		}

	free((void *) pWallData);

	return 0;
	}
#endif /* STANDALONE */

/* wallsio.c */

#ifdef DEM

#include <math.h> /* for log10() */
#include "collision.h" /* includes pkd.h and so dem.h */
#include "linalg.h"

#define MIN_LOG_OVERLAP (-6.0)
#define MAX_LOG_OVERLAP ( 0.0)

#define MIN_COS_A (-1.0)
#define MAX_COS_A ( 1.0)

#define MIN_LOG_S (-5.0)
#define MAX_LOG_S ( 1.0)

#ifndef SQ
double SQ(double);
#define SQ(x) ((x)*(x))
#endif

int demCheckOverlapPoint(const Vector vRelPos,double dRadSq)
{
	return vectorMagSq(vRelPos) <= dRadSq;
}

#if defined(WALLS) && defined(WALLS_REACT)

void pkdDEMWallsReact(PKD pkd,double *dTotalZForceFromParticles) {
	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	*dTotalZForceFromParticles = 0.0;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		*dTotalZForceFromParticles += p->dZForceOnWalls;
		}
	}

#endif /* defined(WALLS) && defined(WALLS_REACT) */

void pkdDEMZeroSmallMotions(PKD pkd,double dAccCritSq,double dDeltaSq)
{
	/*
	** Zeroes accelerations and velocities (and torques per unit mass
	** and spins) for particles with these quantities below a
	** user-specified threshold.  The strategy is that the
	** acceleration quantities must be oppositely directed to the
	** velocity quantities, and both must have magnitudes below the
	** threshold (a critical acceleration, e.g., a small fraction of
	** the ambient gravity) to be zeroed.  Note: at this point p->v is
	** the half-step velocity, but that's ok because if we zero both
	** p->a and p->v, then the full-step velocity will also be zero.
	** The same is true for p->w.  This routine should be called after
	** DEM forces have been fully calculated and before the closing
	** velocity kick.
	*/

	PARTICLE *p;
	double r2,dp,a2,v2;
	int i,n;

	n = pkdLocal(pkd);
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		r2 = RADIUS(p)*RADIUS(p);
		assert(r2 > 0.0);
		dp = p->a[0]*p->v[0] + p->a[1]*p->v[1] + p->a[2]*p->v[2]; /* dv/dt dot v */
		if (dp < 0.0) {
			a2 = p->a[0]*p->a[0] + p->a[1]*p->a[1] + p->a[2]*p->a[2];
			if (a2 > 0.0 && a2 < dAccCritSq) {
				v2 = p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2];
				if (v2 > 0.0 && v2 < dAccCritSq*dDeltaSq) {
					p->v[0] = p->v[1] = p->v[2] = 0.0;
					p->a[0] = p->a[1] = p->a[2] = 0.0;
					}
				}
			}
		dp = p->wDot[0]*p->w[0] + p->wDot[1]*p->w[1] + p->wDot[2]*p->w[2]; /* dw/dt dot w */
		if (dp < 0.0) {
			a2 = p->wDot[0]*p->wDot[0] + p->wDot[1]*p->wDot[1] + p->wDot[2]*p->wDot[2];
			if (a2 > 0.0 && a2 < dAccCritSq*r2) {
				v2 = p->w[0]*p->w[0] + p->w[1]*p->w[1] + p->w[2]*p->w[2];
				if (v2 > 0.0 && v2 < dAccCritSq*dDeltaSq/r2) {
					p->w[0] = p->w[1] = p->w[2] = 0.0;
					p->wDot[0] = p->wDot[1] = p->wDot[2] = 0.0;
					}
				}
			}
		}
	}

void pkdDEMStats(PKD pkd,DEM_STATS *ds) {
	PARTICLE *p;
	double dWidth,d;
	int i,nLocal = pkdLocal(pkd);

	/* find minimum particle separation & maximum relative overlap speed on this processor */

	ds->fDistMin = FLOAT_MAXVAL;
	ds->fSpeedMax = 0.0;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->DEMStats.fDistMin < ds->fDistMin)
			ds->fDistMin = p->DEMStats.fDistMin;
		if (p->DEMStats.fSpeedMax > ds->fSpeedMax)
			ds->fSpeedMax = p->DEMStats.fSpeedMax;
		}

	assert(ds->fDistMin >= 0.0); /* sanity check */
	assert(ds->fSpeedMax >= 0.0);

	/* initialize histogram bins */

	assert(DEM_NUM_OVERLAP_BINS > 0);

	for (i=0;i<DEM_NUM_OVERLAP_BINS;i++)
		ds->pOverlapHist[i] = 0;

	dWidth = (MAX_LOG_OVERLAP - MIN_LOG_OVERLAP)/DEM_NUM_OVERLAP_BINS;

	/* fill overlap histogram with local particle data */

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->DEMStats.fOverlap > 0.0) {
			d = log10(p->DEMStats.fOverlap);
			if (d < MIN_LOG_OVERLAP)
				++ds->pOverlapHist[0];
			else if (d >= MAX_LOG_OVERLAP) /* > MAX only possible if small particle fully inside large particle */
				++ds->pOverlapHist[DEM_NUM_OVERLAP_BINS - 1];
			else
				++ds->pOverlapHist[(int) ((d - MIN_LOG_OVERLAP)/dWidth)];
			}
		}

	/* repeat for cos(alpha) histogram */

	assert(DEM_NUM_COS_A_BINS > 0);

	for (i=0;i<DEM_NUM_COS_A_BINS;i++)
		ds->pCosAHist[i] = 0;

	dWidth = (MAX_COS_A - MIN_COS_A)/DEM_NUM_COS_A_BINS;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		d = p->DEMStats.fCosA;
		if (d >= MIN_COS_A && d <= MAX_COS_A) {
			if (d == MAX_COS_A)
				++ds->pCosAHist[DEM_NUM_COS_A_BINS - 1];
			else
				++ds->pCosAHist[(int) ((d - MIN_COS_A)/dWidth)];
			}
		}

	/* ditto for S vectors */

	assert(DEM_NUM_S_BINS > 0);

	for (i=0;i<DEM_NUM_S_BINS;i++)
		ds->pSHist[i] = 0;

	dWidth = (MAX_LOG_S - MIN_LOG_S)/DEM_NUM_S_BINS;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->DEMStats.fS2 > 0.0) {
			d = 0.5*log10(p->DEMStats.fS2/SQ(RADIUS(p)));
			if (d < MIN_LOG_S)
				++ds->pSHist[0];
			else if (d >= MAX_LOG_S)
				++ds->pSHist[DEM_NUM_S_BINS - 1];
			else
				++ds->pSHist[(int) ((d - MIN_LOG_S)/dWidth)];
			}
		}
	}

#ifdef AGGS

#include "aggs.h"

void pkdDEMAggsSetMass(PKD pkd,int iAggIdx,FLOAT fMass)
{
	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx)
			p->fMassAgg = fMass;
		}
	}

#endif /* AGGS */

#endif /* DEM */

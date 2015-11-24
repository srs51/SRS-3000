/* aggs.c
 * 
 * PKDGRAV Source Code for Aggregate Handling
 * 
 * Author: Kenneth W. Flynn
 *         flynnk@astro.umd.edu
 * Mods:   Derek C. Richardson
 *         dcr@astro.umd.edu
 *
 * Modified: 01/28/01; DCR: 07/10/02, 5/29/03, 7/14/05
 */
 
#ifdef AGGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "aggs.h"
#include "collision.h"

#ifdef AGGS_IN_PATCH
static struct {
  double dStepTime;
  double dEventTime;
  double dOmega;
  } aggs_in_patch_extras; /* global to aggs.c (not ideal, but far simpler) */
#endif /* AGGS_IN_PATCH */

void pkdAggsFind(PKD pkd,int *iMaxIdx)
{
	/*
	 ** Returns largest aggregate index found on local processor.
	 */

	PARTICLE *p;
	int i,n,iAggIdx;

	*iMaxIdx = -1;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (IS_AGG(p)) {
			iAggIdx = AGG_IDX(p);
			if (iAggIdx > *iMaxIdx)
				*iMaxIdx = iAggIdx;
			}
		}
	}

void pkdAggsCountPart(PKD pkd,int iAggIdx,int *nPart)
{
	/*
	 ** Counts number of particles belonging to aggregate (if any).
	 */

	PARTICLE *p;
	int i,n;

	*nPart = 0;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx)
			++(*nPart);
		}
	}

void pkdAggsMerge(PKD pkd,int iOldIdx,int iNewIdx)
{
	/*
	 ** "Merges" two aggregates by assigning aggregate index of each
	 ** particle in old aggregate to index of new aggregate.  A call
	 ** to msrAggsUpdate() is needed to compute new dynamical
	 ** quantities (COM, spin, etc.; this is done in msrAggsMerge()).
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	/* for now, the color of the particles is preserved */

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iOldIdx) {
			AGG_SET_IDX(p,iNewIdx);
			}
		}
	}

void pkdAggsBackDrift(PKD pkd,int iAggIdx,double dt)
{
	/*
	 ** Drifts aggregate particle space positions back interval dt
	 ** (should be to start of step) so that collision prediction
	 ** (which assumes particle positions are at start of step) will
	 ** work properly.  Also sets SMOOTHACTIVE for all particles in
	 ** aggregate to force recomputation of collision circumstances.
	 ** Particle "accelerations" are taken to be second-order terms in
	 ** velocity expression, as computed in pkdAggsSetSpaceVel(),
	 ** which is called in msrAggsBackDrift() beforehand.
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			if (dt > 0.0) { /* (don't backdrift "into the future" ==> overlap problem) */
				p->r[0] -= (p->v[0] + 0.5*p->omega_v[0]*dt)*dt;
				p->r[1] -= (p->v[1] + 0.5*p->omega_v[1]*dt)*dt;
				p->r[2] -= (p->v[2] + 0.5*p->omega_v[2]*dt)*dt;
				}
			TYPESet(p,TYPE_SMOOTHACTIVE);
			}
		}
	}

#ifdef AGGS_IN_PATCH
void pkdAggsGetCOM(PKD pkd,int iAggIdx,const PATCH_PARAMS *PP,Scalar *m,Vector mr,Vector mv/*,Scalar *mPy*/)
#else
void pkdAggsGetCOM(PKD pkd,int iAggIdx,Scalar *m,Vector mr,Vector mv)
#endif
{
	/*
	 ** Computes contribution (moments) of local particles to center-
	 **  of-mass position and velocity of specified aggregate.  Also
	 **  computes total mass of aggregate.
	 ** As of 9/14/09, it not longer computes the mass-weighted dPy in the
	 **  case of AGGS_IN_PATCH.  This is handled whenever an event might change 
	 **  Py: bounce, frag, merge, break via stress, or kick.
	 */

	PARTICLE *p;
	int i,n;

	*m = 0.0;
	vectorZero(mr);
	vectorZero(mv);

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			*m += p->fMass;
#ifdef AGGS_IN_PATCH
			mr[0] += p->fMass*p->r_unwrapped[0];
			mr[1] += p->fMass*p->r_unwrapped[1];
			mr[2] += p->fMass*p->r_unwrapped[2];
			mv[0] += p->fMass*p->v_unwrapped[0];
			mv[1] += p->fMass*p->v_unwrapped[1];
			mv[2] += p->fMass*p->v_unwrapped[2];
#else
			mr[0] += p->fMass*p->r[0];
			mr[1] += p->fMass*p->r[1];
			mr[2] += p->fMass*p->r[2];
			mv[0] += p->fMass*p->v[0];
			mv[1] += p->fMass*p->v[1];
			mv[2] += p->fMass*p->v[2];
#endif
			}
		}
	}

#ifdef AGGS_IN_PATCH

void pkdAggsInPatchGetRef(PKD pkd,int iAggIdx,Scalar *m_max,Vector r_max)
{
	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	*m_max = -1.0; /* this indicates no non-negative masses found */
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx && p->fMass > *m_max) {
			*m_max = p->fMass;
			vectorCopy(p->r,r_max);
			}
		}
	}

void pkdAggsInPatchGetUnwrapped(PKD pkd,int iAggIdx,const Vector r_ref,double dTime,const PATCH_PARAMS *PP)
{
	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {			
			vectorCopy(p->r,p->r_unwrapped);
			vectorCopy(p->v,p->v_unwrapped);

			/* we're assuming a 2-D patch here, with origin at (0,0,0)... */
			if (p->r[0] - r_ref[0] > 0.5*PP->dWidth) {
				p->r_unwrapped[0] -= PP->dWidth;
				p->r_unwrapped[1] += SHEAR(-1,dTime,PP);
				p->v_unwrapped[1] += 1.5*PP->dOrbFreq*PP->dWidth;
				}
			else if (r_ref[0] - p->r[0] > 0.5*PP->dWidth) {
				p->r_unwrapped[0] += PP->dWidth;
				p->r_unwrapped[1] += SHEAR(1,dTime,PP);
				p->v_unwrapped[1] -= 1.5*PP->dOrbFreq*PP->dWidth;
				}
			/* finally, check newly formed unwrapped y-position wrt
			   reference, since shear could place it anywhere */
			if (p->r_unwrapped[1] - r_ref[1] > 0.5*PP->dLength){
				p->r_unwrapped[1] -= PP->dLength;
				}
			else if (r_ref[1] - p->r_unwrapped[1] > 0.5*PP->dLength){
				p->r_unwrapped[1] += PP->dLength;
			}
		}
	}
}

void pkdAggsInPatchOffset(PKD pkd,int iAggIdx,FLOAT dx,FLOAT dy,FLOAT dvy,int bDoUnwrapped)
{
	/*
	** Moves all particles in agg by dx & dy (with y-component of
	** velocity adjusted by dvy).  Called either by
	** msrAggsInPatchCheckCOM() (in response to center of mass
	** boundary crossing) or by msrDoCollisions() (if two aggregates
	** merge).  If bDoUnwrapped is set, the offsets are only applied
	** to "unwrapped" quantities (which is appropriate for a
	** center-of-mass boundary crossing, since p->r and p->v have
	** already been boundary wrapped).  Otherwise, a call to
	** msrAggsInPatchGetUnwrapped() is needed after returning from
	** this function.

	** RP 10/15/09 update: In general, when in-patch ("wrapped")
	** information is altered, it is necessary to wrap dPy as
	** well, in order to remain consistent with in-patch velocity.
	** This operation is essentialy irrelevant at this time, as
	** particles in aggregates recompute their dPy upon liberation
	** from the agg, and do not use their dPy while in the
	** aggregate.  But if we change strategy later, it is good to
	** have this consistency in place.
	*/

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			if (bDoUnwrapped) {
				p->r_unwrapped[0] += dx;
				p->r_unwrapped[1] += dy;
				p->v_unwrapped[1] += dvy;
			        }
			else {
				p->r[0] += dx;
				p->r[1] += dy;
				p->v[1] += dvy;
				p->dPy -= dvy/3.0; /* RP (10/15/09) */
			        }
			}
		}
	}

#endif /* AGGS_IN_PATCH */

void pkdAggsGetAxesAndSpin(PKD pkd,int iAggIdx,const Vector r_com,
						   const Vector v_com,Matrix I,Vector L,
						   Scalar *rad_max_sq,Scalar *volume)
{
	/*
	 ** Computes contribution of local particles to inertia tensor and
	 ** angular momentum vector relative to center of mass of
	 ** specified aggregate.  Particles belonging to the aggregate
	 ** have a copy of their positions relative to the aggregate COM
	 ** stored in p->r_body (in space* coordinates, not body
	 ** coordinates, since the transformation matrix isn't availble
	 ** yet---we're building it now!; to transform later, see
	 ** pkdAggsToBodyAxes()).  Also computes aggregate "radius" as
	 ** maximum particle distance from center of mass (stored as
	 ** square, rad_sq; take square root in calling routine).
	 */

        /* "patch," if in aggs_in_patch */

	PARTICLE *p;
	Vector r,v;
	Scalar m,R;
	double q,r2;
	int i,n;

	I[0][0] = I[0][1] = I[0][2] = I[1][1] =	I[1][2] = I[2][2] = 0.0;
	vectorZero(L);
	*rad_max_sq = *volume = 0.0;

	n = pkdLocal(pkd);
 
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* get pos & vel wrt COM */
#ifdef AGGS_IN_PATCH
			vectorSub(p->r_unwrapped,r_com,r);
			vectorSub(p->v_unwrapped,v_com,v);
#else
			vectorSub(p->r,r_com,r);
			vectorSub(p->v,v_com,v);
#endif
			m = p->fMass;
			R = RADIUS(p); /* particle radius */
			q = AGGS_PARTICLE_INERTIA_PREFACTOR*R*R; /* for convenience */
			/* add inertia tensor contributions */
			/* note caller must fill symmetric elements; Cf. msrAggsGetAxesAndSpin() */
			I[0][0] += m*(q + r[1]*r[1] + r[2]*r[2]); 
			I[0][1] -= m*r[0]*r[1];
			I[0][2] -= m*r[0]*r[2];
			I[1][1] += m*(q + r[0]*r[0] + r[2]*r[2]);
			I[1][2] -= m*r[1]*r[2];
			I[2][2] += m*(q + r[0]*r[0] + r[1]*r[1]);
			/* add angular momentum contributions, L = m (r x v) + I w */
			/* note w for aggregate particles should equal w of aggregate as a whole */
			L[0] += m*(r[1]*v[2] - r[2]*v[1] + p->w[0]*q);
			L[1] += m*(r[2]*v[0] - r[0]*v[2] + p->w[1]*q);
#ifdef AGGS_IN_PATCH
			/* need to account for patch rotation */
			L[2] += m*(r[0]*v[1] - r[1]*v[0] + (p->w[2] - pkd->PP->dOrbFreq)*q); /* RP-DEBUG: hold on; we use patch coords to make spin in patch frame.  Use particle spin - from space frame - and subtract Omega in z to put spin into patch.  BUT do we adjust spin orientation into patch?  (Looks like not!) */
#else
			L[2] += m*(r[0]*v[1] - r[1]*v[0] + p->w[2]*q);
#endif
			/* check to see if aggregate maximum "radius" should be updated */
			r2 = vectorMagSq(r);
			if (r2 > *rad_max_sq)
				*rad_max_sq = r2;
			/* increment "volume" */
			*volume += R*R*R;
			/* store pos wrt COM */
			vectorCopy(r,p->r_body);
			}
		}
	}

void pkdAggsSetBodyPos(PKD pkd,int iAggIdx,Matrix spaceToBody)
{
	/*
	 ** Transforms positions of local particles (belonging
	 ** to specified aggregate) from space to body coordinates
	 ** relative to center of mass.  Note pkdAggsGetAxes() must
	 ** be called first.
	 */

	PARTICLE *p;
	Vector tmp;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			matrixTransform(spaceToBody,p->r_body,tmp);
			vectorCopy(tmp,p->r_body);
			}
		}
	}

#ifdef AGGS_IN_PATCH
void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,const Vector r_com,Matrix lambda,double dTime,const PATCH_PARAMS *PP)
#else
void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,const Vector r_com,Matrix lambda)
#endif
{
	/*
	 ** Transforms positions of local particles (belonging to
	 ** specified aggregate) from body to space coordinates.
	 ** Called by msrAggsAdvance() during drift step.
	 */

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* transform positions from body to space coords */
			matrixTransform(lambda,p->r_body,p->r);
			/* add center of mass component */
			vectorAdd(p->r,r_com,p->r);
#ifdef AGGS_IN_PATCH
			vectorCopy(p->r,p->r_unwrapped);
			/* y-vel adjusted in pkdSetSpaceVel() if needed */
			if (p->r[0] > 0.5*PP->dWidth) {
				p->r[0] -= PP->dWidth;
				p->r[1] += SHEAR(-1,dTime,PP);
				}
			else if (p->r[0] < -0.5*PP->dWidth) {
				p->r[0] += PP->dWidth;
				p->r[1] += SHEAR(1,dTime,PP);
				}
			if (p->r[1] > 0.5*PP->dLength)
				p->r[1] -= PP->dLength;
			else if (p->r[1] < -0.5*PP->dLength)
				p->r[1] += PP->dLength;
#endif /* AGGS_IN_PATCH */
			}
		}
	}

void pkdAggsSetSpaceVel(PKD pkd,int iAggIdx,const Vector v_com,
						const Vector omega,Matrix lambda)
{
	/*
	 ** Computes space velocities of local particles (belonging
	 ** to specified aggregate) to 2nd order.  Called after a
	 ** kick by msrAggsKick() and before back drifting by
	 ** msrAggsBackDrift().
	 */

	PARTICLE *p;
	Vector v,u;
	int i,n;

#ifdef AGGS_IN_PATCH
	PATCH_PARAMS *PP = pkd->PP;
	Vector omega_patch;
#endif

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
#ifdef AGGS_IN_PATCH
			/*
			** For this case, we need to first transform to the space
			** frame, then to the patch frame...
			*/
			matrixTransform(lambda,omega,omega_patch); /* --> space */
			omega_patch[2] -= PP->dOrbFreq; /* --> patch */ 
			                                /* RP-DEBUG: ROTATION 
							   based on time from start
							   of the simulation provided 
							   by transform via lambda. 
							*/
			matrixTransform(lambda,p->r_body,u); /* not safe to use p->r */
			vectorCross(omega_patch,u,v);
			vectorAdd(v,v_com,p->v);
			vectorCopy(p->v,p->v_unwrapped);
			p->v[1] -= 1.5*PP->dOrbFreq*(p->r[0] - p->r_unwrapped[0]);
			vectorCross(omega_patch,v,p->omega_v);
#else
			/* compute aggregate spin component of particle velocities */
			vectorCross(omega,p->r_body,v);
			/* transform to space frame */
			matrixTransform(lambda,v,p->v);
			/* add center of mass component */
			vectorAdd(p->v,v_com,p->v);
			/* compute 2nd-order (centripetal) term */
			vectorCross(omega,v,u);
			/* transform & store in particle structure */
			matrixTransform(lambda,u,p->omega_v);
#endif /* !AGGS_IN_PATCH */
			}
		}
	}

void pkdAggsSetSpaceSpins(PKD pkd,int iAggIdx,const Vector omega)
{
	/*
	** Sets particle spins to aggregate spin.  NOTE: omega is the
	** *space-frame* spin of the aggregate, computed in
	** msrAggsSetSpaceSpins().
	*/

	PARTICLE *p;
	int i,n;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx)
			vectorCopy(omega,p->w);
		}
	}

static void aggsReleaseParticle(PARTICLE *p)
{
	/*
	** Alters particle data to indicate that particle no longer
	** belongs to an aggregate.  Currently this is accomplished by
	** setting the original index to a non-negative integer (for now,
	** INT_MAX) and changing the color back to 3 (green).  In the
	** future, perhaps something more sophisticated will be done here.
	*/

	assert(p != NULL);
	p->iOrgIdx = INT_MAX; /* for now */
	if (ALLOW_COLOR_CHANGE(p))
		p->iColor = 3/*PLANETESIMAL*/; /* ditto */
#ifdef AGGS_IN_PATCH
	/* derotate for length of simulation */
	aggsRotateForPatch(p->w,-1.0);
#endif
	}

#ifdef AGGS_IN_PATCH
void pkdAggsDelete(PKD pkd,int iAggIdx,double dStepTime,double dEventTime,int *bFound)
#else
void pkdAggsDelete(PKD pkd,int iAggIdx,int *bFound)
#endif
{
	PARTICLE *p;
	int i,n;

#ifdef AGGS_IN_PATCH
	/* 
	   Reminder: the info in pkd->PP isn't reliable after a domain
	   decomposition, and must be recopied into PP at the pst level.
	*/
	aggs_in_patch_extras.dStepTime = dStepTime;
	aggs_in_patch_extras.dEventTime = dEventTime;
	aggs_in_patch_extras.dOmega = pkd->PP->dOrbFreq;
#endif

	*bFound = 0;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
#ifdef AGGS_IN_PATCH
		    /* calculate sole particle's Py before 'deleting' */
		    p->dPy = p->v[1] + 2.0*pkd->PP->dOrbFreq*
		      (p->r[0] + p->v[0]*((pkd->PP->dDelta/2.0) - dEventTime)); /*RP-DEBUG-dPy*/
#endif
	            aggsReleaseParticle(p);
		    *bFound = 1;
		    return;
		    }
		}
	}

void pkdAggsGetAccel(PKD pkd,int iAggIdx,Scalar *m,Vector ma)
{
	/*
	 ** Computes contribution (moments) of local particles to center
	 ** of mass acceleration of specified aggregate.  Must be called
	 ** after computing interparticle gravitational accelerations.
	 */

	PARTICLE *p;
	int i,n;

	*m = 0.0;
	vectorZero(ma);

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			*m += p->fMass; /* used as a check in msrAggsGetAccel() */
			ma[0] += p->fMass*p->a[0];
			ma[1] += p->fMass*p->a[1];
			ma[2] += p->fMass*p->a[2];
			}
		}
	}

#ifndef UNIT_CONV
#define UNIT_CONV
#define L_SCALE 1.49597870e11 /* 1 AU in m */
#define M_SCALE 1.9891e30 /* M_Sun in kg */
#define T_SCALE (3.15581497632e7/M_TWO_PI) /* 1 sid yr/2 pi in s */
#endif

void pkdAggsCheckStress(PKD pkd,int iAggIdx,
						Scalar mass,Scalar rad_max,Scalar rad_eff,
						const Vector r_com,const Vector a_com,
						const Vector omega,STRENGTH_PARAMS *SP,
						int *nLost,int *nLeft)
{
	/*
	** Checks each particle in agg to see if acceleration wrt com
	** exceeds threshold, in which case particle is liberated.
	*/

	/* cgs <--> pkdgrav conversion factors */

	const double dLS = 1.49597870e13; /* 1 AU in cm */
	const double dMS = 1.9891e33; /* 1 M_Sun in g */
	const double dTS = 3.15581497632e7/(2.0*M_PI); /* 1 sid yr/2pi in s */

	PARTICLE *p;
	Vector r,r_hat,a,v,a_cen,a_tot,a_tns,a_shr;
	Scalar r_mag;
	double dTensileStrength0,dShearStrength0,dConvFac,dScaleFac,dTensileStrength,dShearStrength;
	double dTensileStress,dShearStress;
	int i,n;

#ifdef AGGS_IN_PATCH
	aggs_in_patch_extras.dStepTime = pkd->dTime;
	aggs_in_patch_extras.dEventTime = pkd->PP->dDelta; /* Revised on 10/13/2009.  Was 0.0 */
	aggs_in_patch_extras.dOmega = pkd->PP->dOrbFreq;
#endif

	/*
	** Compute the aggregate strengths (in dyne/cm^2) based on rough
	** measure of size.  Only need to do this once at the beginning.
	*/

	dTensileStrength0 = (SP->dTensileCoef == DBL_MAX ? DBL_MAX :
						 (SP->dTensileCoef == 0.0 ? 0.0 :
						  SP->dTensileCoef*pow(rad_max*dLS,SP->dTensileExp)));

	dShearStrength0 = (SP->dShearCoef == DBL_MAX ? DBL_MAX :
					   (SP->dShearCoef == 0.0 ? 0.0 :
						SP->dShearCoef*pow(rad_max*dLS,SP->dShearExp)));

	/*
	** Compute maximum allowed average differential acceleration,
	** putting the result back in pkdgrav units.  To be confusing, we
	** continue to call this a "strength", even though the units of
	** strength are force per unit area, not acceleration.
	*/

	dConvFac = dLS*dTS*dTS/dMS; /* from dyne/cm^2 to pkdgrav units */

	/*DEBUG! 7/23/07 -- for now, let's use the particle radius and
	  mass, instead of aggregate effective radius and total mass, to
	  see how that affects the behaviour.  So, we block out the next
	  few lines and put the new calculation in the loop over particles
	  below...*/

#ifdef OLD_STRENGTH_STUFF
	dScaleFac = M_PI*rad_eff*rad_eff/mass; /* to get max acceleration */
	dTensileStrength = dTensileStrength0*dConvFac*dScaleFac;
	dShearStrength = dShearStrength0*dConvFac*dScaleFac;
#endif

	/* now check each particle in turn */

	*nLost = *nLeft = 0;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
			/* particle position relative to com */
#ifdef AGGS_IN_PATCH
			vectorSub(p->r_unwrapped,r_com,r); /* use true distance to COM */
#else
			vectorSub(p->r,r_com,r);
#endif
			/* unit radial vector from com to particle */
			if ((r_mag = vectorMag(r)) == 0.0){
			        ++(*nLeft);
				continue; /* particle located at center of mass! */
			}
			vectorScale(r,1.0/r_mag,r_hat);
			/* differential acceleration wrt com */
			vectorSub(p->a,a_com,a);
#ifdef AGGS_IN_PATCH
			{
				/* add tidal & centrifugal terms to *differential* acceleration (Cf. Dan's Eq. (18)) */
				a[0] += 3.0*pkd->PP->dOrbFreq*pkd->PP->dOrbFreq*r[0];
				a[2] -= pkd->PP->dOrbFreqZ2*r[2];
				}
#endif /* AGGS_IN_PATCH */
			/* now add centrifugal term due to spin around mass center */
			vectorCross(omega,r,v);
			vectorCross(omega,v,a_cen);
			vectorSub(a,a_cen,a_tot); /* sub: centripetal --> centrifugal */
			/* compute stress */
			dTensileStress = vectorDot(a_tot,r_hat);
			vectorScale(r_hat,dTensileStress,a_tns); /* tensile stress vec */
			vectorSub(a_tot,a_tns,a_shr); /* shear stress vec */
			dShearStress = vectorMag(a_shr);
#ifndef OLD_STRENGTH_STUFF /*DEBUG! see 7/23/07 comment above*/
			dScaleFac = M_PI*RADIUS(p)*RADIUS(p)/p->fMass; /* to get max acceleration */
			dTensileStrength = dTensileStrength0*dConvFac*dScaleFac;
			dShearStrength = dShearStrength0*dConvFac*dScaleFac;
#endif
			/* (negative [compressional] stress ignored) */
			if (dTensileStress > dTensileStrength ||
				dShearStress > dShearStrength) {
#ifdef AGGS_IN_PATCH
				/* Remaining agg fragment's Py is treated in msrAggsUpdate() */
				p->dPy = p->v[1] + 2.0*pkd->PP->dOrbFreq*
				  (p->r[0] - p->v[0]*(0.5*pkd->PP->dDelta)); /*RP-DEBUG-dPy*/
#endif /* AGGS_IN_PATCH */
				aggsReleaseParticle(p);
				++(*nLost);
				}
			else
				++(*nLeft);
			}
		}
	}

void pkdAggsGetTorque(PKD pkd,int iAggIdx,const Vector r_com,const Vector a_com,
					  Vector torque)
{
	/*
	 ** Computes contribution of local particles to torque of
	 ** specified aggregate in space coordinates relative to center of
	 ** mass.  Note that the COM position is passed, rather than the
	 ** rotation matrix, since it's simpler this way...
	 */

	PARTICLE *p;
	Vector da,dr,cross,dN;
	int i,n;

	vectorZero(torque);

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (AGG_IDX(p) == iAggIdx) {
#ifdef AGGS_IN_PATCH
			vectorSub(p->r_unwrapped,r_com,dr);/* use true distance to COM */
#else
			vectorSub(p->r,r_com,dr); /* or could xform body to space */
#endif
			vectorSub(p->a,a_com,da);
			vectorCross(dr,da,cross);
			vectorScale(cross,p->fMass,dN);
			/* N += m (r x a) */
			vectorAdd(torque,dN,torque);
			}
		}
	}

#ifdef AGGS_IN_PATCH

void aggsRotateForPatch(Vector w,double dDirection)
{
	/*
	** Rotates vector w (usually spin) through angle Omega * time,
	** where time is the current absolute time since the start of the
	** simulation, to account for patch rotation relative to space
	** frame.  Needed to get single particle spin orientations right
	** for collisions.  Similar to msrAggsInPatchRotateAxes().
	*/

	Matrix mPatchRot;
	Vector v;
	double dTime,alpha,c,s;

	dTime = aggs_in_patch_extras.dStepTime + aggs_in_patch_extras.dEventTime;
	if (dTime == 0.0)
		return;
	alpha = - aggs_in_patch_extras.dOmega*dTime*dDirection;
	c = cos(alpha);
	s = sin(alpha);
	matrixIdentity(mPatchRot);
	mPatchRot[0][0] = c;
	mPatchRot[0][1] = -s;
	mPatchRot[1][0] = s;
	mPatchRot[1][1] = c;
	matrixTransform(mPatchRot,w,v);
	vectorCopy(v,w);
	}

#endif /* AGGS_IN_PATCH */

static void aggsSingleParticleToAgg(const COLLIDER *c,Aggregate *a)
{
	double inertia;

	a->bAssigned = 1;
	a->mass = c->fMass;
	vectorCopy(c->r,a->r_com);
	vectorCopy(c->v,a->v_com);
	vectorZero(a->a_com); /* unused */
	vectorZero(a->torque); /* unused */
	vectorCopy(c->w,a->omega);
#ifdef AGGS_IN_PATCH
	a->dPy_com = c->dPy;
	/* need to orient the "aggregate" spin vector to account for patch rotation */
	aggsRotateForPatch(a->omega,1.0);
#endif /* AGGS_IN_PATCH */
	inertia = AGGS_PARTICLE_INERTIA_PREFACTOR*c->fMass*c->fRadius*c->fRadius;
	vectorSet(a->moments,inertia,inertia,inertia);
	matrixIdentity(a->lambda);
	a->dLastUpdate = 0.0; /* unused */
	}

static void aggsAggToSingleParticle(const Aggregate *a,COLLIDER *c)
{
	/* id, fRadius, iColor, dt, iRung, bTinyStep, agg unchanged */

	c->fMass = a->mass; /* unchanged in aggsBounce() */
	vectorCopy(a->r_com,c->r); /* unchanged in aggsBounce() */
	vectorCopy(a->v_com,c->v);
	vectorCopy(a->omega,c->w);
#ifdef AGGS_IN_PATCH
	c->dPy = a->dPy_com;
	/* need to de-rotate spin vector orientation back to space frame */
	aggsRotateForPatch(c->w,-1.0);
#endif /* AGGS_IN_PATCH */
	}

#ifdef AGGS_IN_PATCH
static void aggsBounce(const COLLIDER *pc1,const COLLIDER *pc2,
					   double dEpsN,COLLIDER *pcOut[],int *pnOut,
		       PATCH_PARAMS *PP, double dCollisionTime)
#else
static void aggsBounce(const COLLIDER *pc1,const COLLIDER *pc2,
					   double dEpsN,COLLIDER *pcOut[],int *pnOut)
#endif
{
	/*
	 ** Collision between aggregates of spheres, from Richardson 1995:
	 ** 	dv1 = gamma (1 + en) (m2/M) un n;
	 ** 	dv2 = - (m1/m2) dv1;
	 ** 	dw1 = m1 I1inv (c1 cross dv1);
	 **		dw2 = - m1 I2inv (c2 cross dv1);
	 ** where M = m1 + m2, un = normal of total relative velocity at
	 ** impact point, and gamma depends on the reduced mass and elements
	 ** of I1inv and I2inv in the ntp basis.
	 **
	 ** NOTE: method assumes NO sliding friction (et = 1).
	 */

	Aggregate *agg1,*agg2;
	Matrix ntpT,ntp,mtmp1,mtmp2,Ibody,spaceToBody1,spaceToBody2,Ispace1,Ispace2,a,b;
	Vector v,n,a1,a2,b1,b2,c1,c2,w1,w2,s1,s2,s,u,t,p,dv1,dv2,dw1,dw2,vtmp;
	double en,m1,m2,M,mu,R1,R2,R,un,c1t,c1p,c2t,c2p,gamma;

	*pnOut = 2;
	*pcOut = (COLLIDER *) malloc((*pnOut)*sizeof(COLLIDER));
	assert(*pcOut != NULL);

	(*pcOut)[0] = *pc1; /* struct copy */
	(*pcOut)[1] = *pc2;

	/*
	 ** Handy pointers.
	 ** NOTE: point to *output* here because we want to change values.
	 */

	agg1 = &((*pcOut)[0].agg);
	agg2 = &((*pcOut)[1].agg);

	/*
	 ** Rather than handle single particles as special cases, we turn
	 ** them into single-particle aggregates for the purpose of solving
	 ** the restitution equations.  They're turned back into single
	 ** particles at the end.
	 */

	if (!COLLIDER_IS_AGG(pc1))
		aggsSingleParticleToAgg(pc1,agg1);
	if (!COLLIDER_IS_AGG(pc2))
		aggsSingleParticleToAgg(pc2,agg2);

	/* sanity check */

	assert(agg1->bAssigned);
	assert(agg2->bAssigned);

	/* convenient shorthand */

	en = dEpsN;
	m1 = agg1->mass;
	m2 = agg2->mass;
	M = m1 + m2;
	mu = m1*m2/M;
	R1 = pc1->fRadius;
	R2 = pc2->fRadius;
	R = R1 + R2;

	/* v: relative linear velocity of agg centers of mass */

	vectorSub(agg2->v_com,agg1->v_com,v);

#ifdef AGGS_IN_PATCH
	/* RP (10-22-09): Check if both aggregate COM's are on the same side of 
	    the patch; adjust velocity, if needed.
	   'While' is necessary due to the rare occassions when the aggregate
	    COMs are *2* patch widths apart.
	*/
	{
	  Vector r;
	  vectorSub(agg2->r_com,agg1->r_com,r);
	  while(r[0] > 0.5*PP->dWidth) 
	    { 
	      v[1] += 1.5*PP->dOrbFreq*PP->dWidth; 
	      r[0] -= PP->dWidth;
	    }
	  while(-r[0] > 0.5*PP->dWidth) { 
	    v[1] -= 1.5*PP->dOrbFreq*PP->dWidth; 
	    r[0] += PP->dWidth;
	  }
	}/*RP-DEBUG-dPy*/
#endif

	/*
	 ** n: vector perpendicular to tangent plane at impact site, pointing
	 ** from contact sphere of aggregate 1 to contact sphere of aggregate 2.
	 */

	vectorSub(pc2->r,pc1->r,n);

#ifdef AGGS_IN_PATCH
	/*RP-DEBUG Colliding particles may be on opposite sides of the patch; 
	  now we must account for patch boundary conditions for direction of normal vector */
	/* NOTE that velocity adjustments between COMs are handled above, and single particles below */
	if (n[0] > 0.5*PP->dWidth) {
	  n[0] -= PP->dWidth;
	  n[1] += SHEAR(-1,dCollisionTime,PP);
	}
	else if (n[0] < -0.5*PP->dWidth) {
	  n[0] += PP->dWidth;
	  n[1] += SHEAR(1,dCollisionTime,PP);
	}
	if (n[1] > 0.5*PP->dLength) {
	  n[1] -= PP->dLength;
	}
	else if (n[1] < -0.5*PP->dLength) {
	  n[1] += PP->dLength;
	}

#endif /* AGGS_IN_PATCH */

	/*
	 ** a1,a2: positions of colliding spheres relative to agg centers of mass.
	 ** Note agg rotation during drift will cause some error here...
	 */

	vectorSub(pc1->r,agg1->r_com,a1);
	vectorSub(pc2->r,agg2->r_com,a2);

#ifdef AGGS_IN_PATCH
	/*RP-DEBUG Now we must again account for patch boundary conditions */
	/* Note that velocity is adjusted above */
	if (a1[0] > 0.5*PP->dWidth) {
	  a1[0] -= PP->dWidth;
	  a1[1] += SHEAR(-1,dCollisionTime,PP);
	}
	else if (a1[0] < -0.5*PP->dWidth) {
	  a1[0] += PP->dWidth;
	  a1[1] += SHEAR(1,dCollisionTime,PP);
	}
	if (a1[1] > 0.5*PP->dLength) {
	  a1[1] -= PP->dLength;
	}
	else if (a1[1] < -0.5*PP->dLength) {
	  a1[1] += PP->dLength;
	}

	if (a2[0] > 0.5*PP->dWidth) {
	  a2[0] -= PP->dWidth;
	  a2[1] += SHEAR(-1,dCollisionTime,PP);
	}
	else if (a2[0] < -0.5*PP->dWidth) {
	  a2[0] += PP->dWidth;
	  a2[1] += SHEAR(1,dCollisionTime,PP);
	}
	if (a2[1] > 0.5*PP->dLength) {
	  a2[1] -= PP->dLength;
	}
	else if (a2[1] < -0.5*PP->dLength) {
	  a2[1] += PP->dLength;
	}

#endif /* AGGS_IN_PATCH */

	/* b1,b2: position of impact site relative to colliding sphere centers */

	vectorScale(n, R1/R,b1);
	vectorScale(n,-R2/R,b2);

	/* c1,c2: position vectors of impact site relative to agg centers of mass */

	vectorAdd(a1,b1,c1);
	vectorAdd(a2,b2,c2);

	/* w1,w2: angular velocities of aggs in space frame */

	matrixTransform(agg1->lambda,agg1->omega,w1);
	matrixTransform(agg2->lambda,agg2->omega,w2);

#ifdef AGGS_IN_PATCH
	{
		/*
		** The spins stored in the aggregate structures are body-frame
		** spins.  The orientation transform matrix converts these to
		** space-frame spins.  We then need to subtract Omega z-hat to
		** finally convert to patch-frame spins.  We will need to
		** reverse this process at the end.  The spins of single
		** (no-agg) particles are stored in the space frame, but since
		** we don't store their orientations, we also need to rotate
		** their spin axes (and derotate them later)---this is done
		** via the calls to aggsSingleParticleToAgg() and
		** aggsAggToSingleParticle() above.  Note the transformation
		** matrix for single particles is just the identity matrix.
		*/

		w1[2] -= aggs_in_patch_extras.dOmega;
		w2[2] -= aggs_in_patch_extras.dOmega;
		}
#endif /* AGGS_IN_PATCH */

	/* s1,s2: spin velocity of agg at impact site */

	vectorCross(w1,c1,s1);
	vectorCross(w2,c2,s2);

	/* s: relative spin velocity at impact site */

	vectorSub(s2,s1,s);

	/* u: total relative velocity at impact site */

	vectorAdd(v,s,u);

	/* construct ntp basis */

	vectorGetBasis(n,t,p);

	/* un: normal component of u */

	un = vectorDot(u,n);

	/* c1t,c1p,c2t,c2p: transverse components of c1,c2 */

	c1t = vectorDot(c1,t);
	c1p = vectorDot(c1,p);
	c2t = vectorDot(c2,t);
	c2p = vectorDot(c2,p);

	/* get inverse inertia tensors wrt ntp basis */

	vectorCopy(n,ntpT[0]); /* ntpT: matrix whose rows are n, t, and p */
	vectorCopy(t,ntpT[1]);
	vectorCopy(p,ntpT[2]);

	matrixTranspose(ntpT,ntp); /* ntp: matrix whose columns are n, t, and p */

	matrixDiagonal(agg1->moments,Ibody); /* inertia tensor in body frame */
	matrixMultiply(agg1->lambda,Ibody,mtmp1);
	matrixTranspose(agg1->lambda,spaceToBody1);
	matrixMultiply(mtmp1,spaceToBody1,Ispace1); /* Ispace1: inertia tensor of agg 1 in space frame */
	matrixMultiply(ntpT,Ispace1,mtmp1);
	matrixMultiply(mtmp1,ntp,mtmp2); /* inertia tensor in ntp basis */
	matrixInverse(mtmp2,a); /* a: inverse of inertia tensor of agg 1 in ntp basis */

	matrixDiagonal(agg2->moments,Ibody); /* inertia tensor in body frame */
	matrixMultiply(agg2->lambda,Ibody,mtmp1);
	matrixTranspose(agg2->lambda,spaceToBody2);
	matrixMultiply(mtmp1,spaceToBody2,Ispace2); /* Ispace2: inertia tensor of agg 2 in space frame */
	matrixMultiply(ntpT,Ispace2,mtmp1);
	matrixMultiply(mtmp1,ntp,mtmp2); /* inertia tensor in ntp basis */
	matrixInverse(mtmp2,b); /* b: inverse of inertia tensor of agg 2 in ntp basis */

	/* useful factor */

	gamma = 1.0/(1.0 + mu*(a[1][1]*c1p*c1p - 2.0*a[1][2]*c1p*c1t + a[2][2]*c1t*c1t +
						   b[1][1]*c2p*c2p - 2.0*b[1][2]*c2p*c2t + b[2][2]*c2t*c2t));

	/* verbose agg conservation check...
	{
		Vector P1,P2,P,L1t,L1rb,L1r,L1,L2t,L2rb,L2r,L2,L;
		Matrix inertia;

		printf("r1 = (%g,%g,%g) v1 = (%g,%g,%g) w1 = (%g,%g,%g)\n",
			   agg1->r_com[0],agg1->r_com[1],agg1->r_com[2],
			   agg1->v_com[0],agg1->v_com[1],agg1->v_com[2],
			   agg1->omega[0],agg1->omega[1],agg1->omega[2]);
		printf("r2 = (%g,%g,%g) v2 = (%g,%g,%g) w2 = (%g,%g,%g)\n",
			   agg2->r_com[0],agg2->r_com[1],agg2->r_com[2],
			   agg2->v_com[0],agg2->v_com[1],agg2->v_com[2],
			   agg2->omega[0],agg2->omega[1],agg2->omega[2]);
		vectorScale(agg1->v_com,agg1->mass,P1);
		vectorScale(agg2->v_com,agg2->mass,P2);
		vectorAdd(P1,P2,P);
		printf("lin mom before = %g %g %g\n",P[0],P[1],P[2]);
		vectorCross(agg1->r_com,agg1->v_com,L1t);
		vectorScale(L1t,agg1->mass,L1t);
		matrixDiagonal(agg1->moments,inertia);
		matrixTransform(inertia,agg1->omega,L1rb);
		matrixTransform(agg1->lambda,L1rb,L1r);
		printf("ang mom before L1t = %g %g %g\n",L1t[0],L1t[1],L1t[2]);
		printf("ang mom before L1r = %g %g %g\n",L1r[0],L1r[1],L1r[2]);
		vectorAdd(L1t,L1r,L1);
		vectorCross(agg2->r_com,agg2->v_com,L2t);
		vectorScale(L2t,agg2->mass,L2t);
		matrixDiagonal(agg2->moments,inertia);
		matrixTransform(inertia,agg2->omega,L2rb);
		matrixTransform(agg2->lambda,L2rb,L2r);
		printf("ang mom before L2t = %g %g %g\n",L2t[0],L2t[1],L2t[2]);
		printf("ang mom before L2r = %g %g %g\n",L2r[0],L2r[1],L2r[2]);
		vectorAdd(L2t,L2r,L2);
		vectorAdd(L1,L2,L);
		printf("ang mom before = %g %g %g\n",L[0],L[1],L[2]);
		}
	*/

	/* compute final velocities and spins */

	vectorScale(n,gamma*(1 + en)*(m2/M)*un,dv1);
	vectorAdd(agg1->v_com,dv1,agg1->v_com);

	vectorScale(dv1,-m1/m2,dv2);
	vectorAdd(agg2->v_com,dv2,agg2->v_com);

#ifdef AGGS_IN_PATCH
	/* remove patch spin */
	w1[2] += aggs_in_patch_extras.dOmega;
	w2[2] += aggs_in_patch_extras.dOmega;
#endif

	vectorCross(c1,dv1,vtmp);
	matrixInverse(Ispace1,a); /* no longer in ntp basis */
	matrixTransform(a,vtmp,dw1);
	vectorScale(dw1,m1,dw1);
	vectorAdd(w1,dw1,w1);
	matrixTransform(spaceToBody1,w1,agg1->omega); /* new spin in body frame */

	vectorCross(c2,dv1,vtmp);
	matrixInverse(Ispace2,b); /* no longer in ntp basis */
	matrixTransform(b,vtmp,dw2);
	vectorScale(dw2,-m1,dw2);
	vectorAdd(w2,dw2,w2);
	matrixTransform(spaceToBody2,w2,agg2->omega); /* new spin in body frame */

	/* verbose agg conservation check...
	{
		Vector P1,P2,P,L1t,L1rb,L1r,L1,L2t,L2rb,L2r,L2,L;
		Matrix inertia;

		printf("r1 = (%g,%g,%g) v1 = (%g,%g,%g) w1 = (%g,%g,%g)\n",
			   agg1->r_com[0],agg1->r_com[1],agg1->r_com[2],
			   agg1->v_com[0],agg1->v_com[1],agg1->v_com[2],
			   agg1->omega[0],agg1->omega[1],agg1->omega[2]);
		printf("r2 = (%g,%g,%g) v2 = (%g,%g,%g) w2 = (%g,%g,%g)\n",
			   agg2->r_com[0],agg2->r_com[1],agg2->r_com[2],
			   agg2->v_com[0],agg2->v_com[1],agg2->v_com[2],
			   agg2->omega[0],agg2->omega[1],agg2->omega[2]);
		vectorScale(agg1->v_com,agg1->mass,P1);
		vectorScale(agg2->v_com,agg2->mass,P2);
		vectorAdd(P1,P2,P);
		printf("lin mom  after = %g %g %g\n",P[0],P[1],P[2]);
		vectorCross(agg1->r_com,agg1->v_com,L1t);
		vectorScale(L1t,agg1->mass,L1t);
		matrixDiagonal(agg1->moments,inertia);
		matrixTransform(inertia,agg1->omega,L1rb);
		matrixTransform(agg1->lambda,L1rb,L1r);
		printf("ang mom  after L1t = %g %g %g\n",L1t[0],L1t[1],L1t[2]);
		printf("ang mom  after L1r = %g %g %g\n",L1r[0],L1r[1],L1r[2]);
		vectorAdd(L1t,L1r,L1);
		vectorCross(agg2->r_com,agg2->v_com,L2t);
		vectorScale(L2t,agg2->mass,L2t);
		matrixDiagonal(agg2->moments,inertia);
		matrixTransform(inertia,agg2->omega,L2rb);
		matrixTransform(agg2->lambda,L2rb,L2r);
		printf("ang mom  after L2t = %g %g %g\n",L2t[0],L2t[1],L2t[2]);
		printf("ang mom  after L2r = %g %g %g\n",L2r[0],L2r[1],L2r[2]);
		vectorAdd(L2t,L2r,L2);
		vectorAdd(L1,L2,L);
		printf("ang mom  after = %g %g %g\n",L[0],L[1],L[2]);
		}
	*/

	/* revert back to single particles as needed */

	if (!COLLIDER_IS_AGG(pc1))
		aggsAggToSingleParticle(agg1,&((*pcOut)[0]));
	if (!COLLIDER_IS_AGG(pc2))
		aggsAggToSingleParticle(agg2,&((*pcOut)[1]));
	}

static void aggsPutColliderInfo(const COLLIDER *c,PARTICLE *p,int iAggIdx)
{
	/* used for merging particles with aggs in pkdAggsDoCollision() */

	vectorCopy(c->r,p->r); /* position at contact */
#ifdef AGGS_IN_PATCH
	vectorCopy(c->v,p->v); /* RP-DEBUG velocity at contact (in case of a ghosted collider */
	p->dPy = c->dPy; /* RP-DEBUG canonical momentum at contact (in case of a ghosted collider */
#endif

	AGG_SET_IDX(p,iAggIdx); /* particle now belongs to this agg */
	if (ALLOW_COLOR_CHANGE(p))
		p->iColor = 4 + iAggIdx%10; /* a quick & dirty way to color aggs... */
	}

#ifdef AGGS_IN_PATCH
void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *pc1,
						const COLLIDER *pc2,int bPeriodic,
						const COLLISION_PARAMS *CP,int iAggNewIdx,
						int *piOutcome,double *dT,
						COLLIDER *cOut,int *pnOut,FLOAT *dx,FLOAT *dy,FLOAT *dvy)
#else
void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *pc1,
						const COLLIDER *pc2,int bPeriodic,
						const COLLISION_PARAMS *CP,int iAggNewIdx,
						int *piOutcome,double *dT,
						COLLIDER *cOut,int *pnOut)
#endif
{
	COLLIDER c1,c2;
	double v2,ve2,dMergeLimit2,dFragLimit2;
	int bReturnOutput,iOutcome,iOverlap,k;

#ifdef AGGS_IN_PATCH
	aggs_in_patch_extras.dStepTime = pkd->dTime;
	aggs_in_patch_extras.dEventTime = dt;
	aggs_in_patch_extras.dOmega = pkd->PP->dOrbFreq;
	*dx = *dy = *dvy = 0.0; /* adjustments for ghost collider */
	double dHalfDelta = 0.5*pkd->PP->dDelta;
#endif

	/* get local copies of collider data for manipulation */

	c1 = *pc1; /* struct copy */
	c2 = *pc2;

	/* verbose collision output...
	printf("COLLISION %i (%i) & %i (%i) (dt = %.16e)\n",
	       pc1->id.iOrder,pc1->id.iOrgIdx,pc2->id.iOrder,pc2->id.iOrgIdx,dt);
	*/

	/* Must be on one of two processors */
	assert(c1.id.iPid == pkd->idSelf || c2.id.iPid == pkd->idSelf); 

	assert(c1.id.iPid == pkd->idSelf || c2.id.iPid == pkd->idSelf); /* otherwise there's a problem!! */

	/* following needed for handling overlap merge condition properly */

	iOutcome = MISS;

	/*
	 ** To prevent overwriting data in parallel, only store results
	 ** in output variables if collider 1 is local to this processor.
	 */

	bReturnOutput = (c1.id.iPid == pkd->idSelf);

	if (bReturnOutput && dT != NULL)
		*dT = 0.0; /* not used (change in energy not computed for aggs) */

#ifdef WALLS
	if (c2.id.iOrder < 0) { /* wall collision */
		assert(bReturnOutput); /* wall is not a particle, so only 1 CPU involved */
		assert(!COLLIDER_STUCK(&c1));
		assert(0); /* not implemented yet */
		/*aggsWallsDoCollision(pkd,CP,&c1,&c2,dt,piOutcome);*/
		assert(*piOutcome == BOUNCE); /* only option supported for now */
		cOut[0] = c1;
		cOut[1] = c2;
		*pnOut = 2;
		return;
		}
#endif /* WALLS */

	/* apply boundary conditions if applicable */

	if (bPeriodic) {
#ifdef AGGS_IN_PATCH
		FLOAT fOffset[3];
		*dvy = pkdApplyBCs(pkd,pkd->dTime/* + dt*/,&c1,&c2,fOffset); /* use dTime + dt here because the aggs have already been advanced to the point of contact *//*DEBUG! only true for aggs, not single particles! not sure how to fix this (it's a matter of consistency: the collision search does not account for future shear; single particles normally adjust for BC shear before advancing to the point of contact; but here aggs are pre-advanced*/
		*dx = fOffset[0];
		*dy = fOffset[1];
		if (COLLIDER_IS_AGG(&c2)) {
			/* need to adjust center-of-mass position and velocity as well */
			c2.agg.r_com[0] += *dx;
			c2.agg.r_com[1] += *dy;
			c2.agg.v_com[1] += *dvy;
			c2.agg.dPy_com -= *dvy/3.0; /* RP-DEBUG: for consistency with pkdApplyBCs() */
			}
#else
		assert(0); /* periodic BCs for AGGS only supported for SLIDING_PATCH */
#endif
		}

	/* handle overlap cases if applicable */

	iOverlap = OverlapIgnore;

	if (dt <= 0.0) {
		iOverlap = CP->iOverlapOption;
		switch (iOverlap) {
		case OverlapIgnore:
		case OverlapBackstep:
			break; /* do nothing */
		case OverlapIsError:
			assert(0);
		case OverlapAdjPos:
			dt = 0.0; /* effect is instantaneous */
			break;
		case OverlapRepel:
			assert(0); /* shouldn't be here */
		case OverlapMerge:
			dt = 0.0;
			break;
		default:
			assert(0);
			}
		}

	/*
	 ** Advance coordinates of non-aggregate particles to impact time.
	 ** (Aggregate particles already advanced in msrAggsAdvance().)
	 ** Note that the aggregate advance step actually integrates the
	 ** Euler equations of motion to the impact time, taking into
	 ** account gravitational torques on each aggregate, whereas
	 ** collision prediction in CheckForCollision() uses a simpler
	 ** expression good to second order assuming the aggregate spin
	 ** vector(s) remain unchanged over the interval.  This means the
	 ** collision circumstances may be slightly off here (particles
	 ** either overlapping or not touching).  To minimize these
	 ** problems, the timestep should be SHORT.  Even so, it may be
	 ** necessary to use one of the overlap options, particularly if
	 ** bouncing is allowed, to circumvent overlap errors.
	 */

	if (!COLLIDER_IS_AGG(&c1)) {
		for (k=0;k<3;k++)
			c1.r[k] += c1.v[k]*dt;
		}
	if (!COLLIDER_IS_AGG(&c2)) {
		for (k=0;k<3;k++)
			c2.r[k] += c2.v[k]*dt;
		}

	/*
	** Compute collision circumstances (relative speed, mutual escape
	** speed at impact), treating aggregates as giant spherical
	** particles for simplicity.
	*/

	/* handle overlap condition, if needed */

	switch (iOverlap) {
	case OverlapIgnore:
	case OverlapBackstep:
		break;
	case OverlapAdjPos:
	  {
		  /* move particles apart along line of centers until just touching */
		  double r_hat[3],r_mag,dr,m_tot;
		  int k;
		  r_mag = 0.0;
		  for (k=0;k<3;k++) {
			  r_hat[k] = c1.r[k] - c2.r[k]; /* points from c2 to c1 */
			  r_mag += r_hat[k]*r_hat[k];
			  }
		  r_mag = sqrt(r_mag); /* separation */
		  assert(r_mag > 0.0); /* can't handle exactly colocated bodies */
		  dr = c1.fRadius + c2.fRadius;
		  assert(r_mag <= dr); /* this is the whole point: they're overlapping */
		  dr -= r_mag; /* expansion factor */
		  assert(dr >= 0.0);
		  m_tot = c1.fMass + c2.fMass;
		  assert(m_tot > 0.0);
		  for (k=0;k<3;k++) {
			  r_hat[k] /= r_mag;
			  c1.r[k] += (c2.fMass/m_tot)*dr*r_hat[k];
			  c2.r[k] -= (c1.fMass/m_tot)*dr*r_hat[k];
			  }
		  /* store new position(s) of aggregate particle(s) */
#ifdef AGGS_IN_PATCH
		  /*DEBUG! we will ignore how this affects dPy for now!*/
#endif
		  if (COLLIDER_IS_AGG(&c1) && c1.id.iPid == pkd->idSelf)
			  pkdPutColliderInfo(pkd,&c1,c2.id.iOrder,&pkd->pStore[c1.id.iIndex],0.0);
		  if (COLLIDER_IS_AGG(&c2) && c2.id.iPid == pkd->idSelf)
			  pkdPutColliderInfo(pkd,&c2,c1.id.iOrder,&pkd->pStore[c2.id.iIndex],0.0);
		  /* handle as normal collision (or miss) from here on */
		  break;
		  }
	case OverlapMerge:
		/* force merge */
		assert(CP->iOutcomes & MERGE);
		iOutcome = MERGE;
		goto merge;
	default:
		assert(0);
		}

	v2 = ve2 = 0.0;

	{ /* RP-DEBUG: only need mutual esc speed if CP->dFragLimit || CP->dMergeLimit > 0! */
		Vector r1,r2;
		double m1,m2,d2=0.0;

		if (COLLIDER_IS_AGG(&c1)) {
			m1 = c1.agg.mass;
			vectorCopy(c1.agg.r_com,r1);
			}
		else {
			m1 = c1.fMass;
			vectorCopy(c1.r,r1);
			}

		if (COLLIDER_IS_AGG(&c2)) {
			m2 = c2.agg.mass;
			vectorCopy(c2.agg.r_com,r2);
			}
		else {
			m2 = c2.fMass;
			vectorCopy(c2.r,r2);
			}

#ifdef AGGS_IN_PATCH
		/* RP-DEBUG: To correctly calculate the mutual escape
		   speed, apply BCs to r1 to account for the
		   possibility that agg COMs are not near the
		   colliders */

		/* The following procedure is identical to
		   pkdAggsInPatchGetUnwrapped(), with r2 providing
		   the "reference" point */

		/* RP-DEBUG:
		    Due to the fact that it is the colliding particles that are 
		    guaranteed to be adjacent -- not the aggregate COMs -- it is
		    possible that when 2 aggregates collide, the aggregate COMs will 
		    be 2 patch widths apart.  If only one correction is made here, 
		    in rare cases (<1% of agg-agg collisions), r1 will be off 
		    by 1 patch width.
		*/

		/* we're assuming a 2-D patch here, with origin at (0,0,0)... */
		if (r1[0] - r2[0] > 0.5*pkd->PP->dWidth) {
		    r1[0] -= pkd->PP->dWidth;
		    r1[1] += SHEAR(-1,pkd->dTime,pkd->PP);

		    /* RP-DEBUG: fixes the double-patch width offset issue... */
		    if (r1[0] - r2[0] > 0.5*pkd->PP->dWidth) {
		      r1[0] -= pkd->PP->dWidth;
		      r1[1] += SHEAR(-1,pkd->dTime,pkd->PP);
		    }
		}
		else if (r2[0] - r1[0] > 0.5*pkd->PP->dWidth) {
		    r1[0] += pkd->PP->dWidth;
		    r1[1] += SHEAR(1,pkd->dTime,pkd->PP);
		    
		    /* Ditto... */
		    if (r2[0] - r1[0] > 0.5*pkd->PP->dWidth) {
		      r1[0] += pkd->PP->dWidth;
		      r1[1] += SHEAR(1,pkd->dTime,pkd->PP);
		    }
		}
		/* finally, check newly formed unwrapped y-position wrt
		   reference, since shear could place it anywhere */
		if (r1[1] - r2[1] > 0.5*pkd->PP->dLength){
		    r1[1] -= pkd->PP->dLength;
		}
		else if (r2[1] - r1[1] > 0.5*pkd->PP->dLength){
		    r1[1] += pkd->PP->dLength;
		}
		
#endif /* AGGS_IN_PATCH */

		for (k=0;k<3;k++) {
			d2 += (r2[k] - r1[k])*(r2[k] - r1[k]);
			/*
			** For aggregates, even though "better" velocities were
			** computed during the advance step immediately prior to
			** processing this collision, we still rely on the
			** velocities used for the original collision prediction
			** for the following computation of the square relative
			** velocity.  Otherwise we would require another step at
			** the master level to replace these collider velocities
			** with the updated ones.  Besides, the collision
			** prediction itself uses those old velocities anyway.
			** For non-aggregates, this is not an issue: the
			** velocities in the collider structs are correct.
			*/
			v2 += (c2.v[k] - c1.v[k])*(c2.v[k] - c1.v[k]);
			}

		/*
		** Since aggregates can have arbitrarily bizarre shapes, it's
		** possible (though unlikely) for particles and/or aggregates
		** to have exactly overlapping mass centers.  So we impose a
		** minimum separation for the escape speed calculation of the
		** sum of the touching particle radii.
		*/

		if (d2 < (c1.fRadius + c2.fRadius)*(c1.fRadius + c2.fRadius))
			d2 = (c1.fRadius + c2.fRadius)*(c1.fRadius + c2.fRadius);

		assert(d2 > 0.0);

		ve2 = 2.0*(m1 + m2)/sqrt(d2);
		}

	/* apply scaling to merge and/or frag limits, if any */

	dMergeLimit2 = CP->dMergeLimit*CP->dMergeLimit;
	if (CP->dMergeLimit > 0.0)
		dMergeLimit2 *= ve2;

	dFragLimit2 = CP->dFragLimit*CP->dFragLimit;
	if (CP->dFragLimit > 0.0)
		dFragLimit2 *= ve2;

	/* determine collision outcome */

	if ((CP->iOutcomes & FRAG) &&
		(COLLIDER_IS_AGG(&c1) || COLLIDER_IS_AGG(&c2)) &&
		v2 > dFragLimit2) {

		/*
		** Liberate each particle from their respective agg (at least
		** 1 collider must be an agg) and allow the particles to
		** bounce (so FRAG also implies BOUNCE in this context).  The
		** master will update what's left of the agg(s).  Note ve2 is
		** not used, and it is assumed that the critical speed for
		** fragmentation exceeds the critical speed for merger (if
		** MERGE is allowed).
		*/

		if (bReturnOutput)
			*piOutcome = FRAG;

		/* release local particles */

		if (COLLIDER_IS_AGG(&c1)) {
#ifdef AGGS_IN_PATCH
		        /* Recompute y-momentum of lost particle to account for fragmentation. */
			c1.dPy = c1.v[1] + 2.0*pkd->PP->dOrbFreq*
			  (c1.r[0] + c1.v[0]*(dHalfDelta - dt)); /*RP-DEBUG-dPy*/

			/*   Remaining aggregate fragment's Py is handled in msrAggsUpdate().
			**     (Note that this happens at msr, after the BOUNCE that is about to
			**     happen to the liberated particle(s).)
			
			     Note -- RP-DEBUG-dPy revision 9/18/09: we edit c1&c2. c[?] are ignored:
			       Once this block is processed, flow goes to 'bounce.'  That code takes c1&c2 as 
			       inputs and creates a copy within aggsBounce(). 
			
			**  Once the bounce is computed, Py computed here will be overwritten, as particle has
			**  new position and velocity info.  (Py is still updated here for consistency, in
			**  case Py is ever used in bounce calculations.)
			*/	
			
#endif /* AGGS_IN_PATCH */
			c1.id.iOrgIdx = INT_MAX; /* indicates collider no longer agg--aggsReleaseParticle() changes pkd struct */
			if (c1.id.iPid == pkd->idSelf) 
				aggsReleaseParticle(&pkd->pStore[c1.id.iIndex]);
		}
		
		if (COLLIDER_IS_AGG(&c2)) {
#ifdef AGGS_IN_PATCH			
			/* recompute y-momentum of lost particle to account for fragmentation */
			c2.dPy = c2.v[1] + 2.0*pkd->PP->dOrbFreq*
			  (c2.r[0] + c2.v[0]*(dHalfDelta - dt)); /*RP-DEBUG-dPy*/
#endif
			c2.id.iOrgIdx = INT_MAX;
			if (c2.id.iPid == pkd->idSelf)
				aggsReleaseParticle(&pkd->pStore[c2.id.iIndex]);
		}
		
		/* now treat particles like bounce */
		
		goto bounce;
		}

	/*
	** If both merging and bouncing are allowed, determine outcome
	** based on rough estimate of mutual escape speed of colliders.
	*/

	merge:

	if (iOutcome == MERGE ||
		(CP->iOutcomes == MERGE ||
		 ((CP->iOutcomes & MERGE) && v2 <= dMergeLimit2))) {

		/*
		 ** Most of the work of merging aggregates is actually done
		 ** at the master level.  Here we're just concerned with
		 ** updating any unaggregated particles to reflect their new
		 ** aggregate member status.

		 RP: note that aggregate Py is also updated in msrAggsUpdate().
		*/

		if (bReturnOutput) {
			*piOutcome = MERGE;
			/* nothing stored in cOut -- master will take care of this */
			*pnOut = 1; /* but pnOut used by pst to show collision found */
			}

		if (COLLIDER_IS_AGG(&c1) && COLLIDER_IS_AGG(&c2)) {
			/* do nothing -- handled in msrAggsMerge() */
		  assert(COLLIDER_AGG_IDX(&c1) != COLLIDER_AGG_IDX(&c2));
		}
		else if (COLLIDER_IS_AGG(&c1) && !COLLIDER_IS_AGG(&c2)) {
			/* add single particle at current position to aggregate */
			if (c2.id.iPid == pkd->idSelf) /* only if particle is local */
				aggsPutColliderInfo(&c2,&pkd->pStore[c2.id.iIndex],COLLIDER_AGG_IDX(&c1));
#ifdef AGGS_IN_PATCH
			/* no further adjustments needed */
			*dx = *dy = *dvy = 0.0;
#endif
			}
		else if (COLLIDER_IS_AGG(&c2) && !COLLIDER_IS_AGG(&c1)) {
			/* ditto */
#ifdef AGGS_IN_PATCH
			/* first move collider 1 so it is adjacent to *unadjusted* collider 2 position */
			c1.r[0] -= *dx;
			c1.r[1] -= *dy;
			c1.v[1] -= *dvy;
			c1.dPy += *dvy/3.0;
			*dx = *dy = *dvy = 0.0;
#endif
			if (c1.id.iPid == pkd->idSelf)
				aggsPutColliderInfo(&c1,&pkd->pStore[c1.id.iIndex],COLLIDER_AGG_IDX(&c2));
			}
		else { /* i.e., !COLLIDER_IS_AGG(&c1) && !COLLIDER_IS_AGG(&c2) */

			/* make new aggregate from single particles at current positions */
			if (c1.id.iPid == pkd->idSelf)
				aggsPutColliderInfo(&c1,&pkd->pStore[c1.id.iIndex],iAggNewIdx);
			if (c2.id.iPid == pkd->idSelf)
				aggsPutColliderInfo(&c2,&pkd->pStore[c2.id.iIndex],iAggNewIdx);
#ifdef AGGS_IN_PATCH
			*dx = *dy = *dvy = 0.0;
#endif
			}
		}
	else if (CP->iOutcomes & BOUNCE) {

		/* bounce */

		if (bReturnOutput)
			*piOutcome = BOUNCE;

	bounce:
		{
			COLLIDER *c;
			double dEpsN = 1.0;
			int i,n;

			/* apply dEpsN rules (see pkdDoCollision()) */
			if (c1.bTinyStep || c2.bTinyStep)
				dEpsN = CP->dCollapseEpsN;
			else if ((CP->iSlideOption == EscVel && v2 < CP->dSlideLimit2*ve2) ||
					 (CP->iSlideOption == MaxTrv && v2 < CP->dSlideLimit2))
				dEpsN = CP->dSlideEpsN;
			else if (v2 > CP->dCrushLimit)
				dEpsN = CP->dCrushEpsN;
			else if (CP->iEpsNOption == ConstEps)
				dEpsN = CP->dEpsN;
			/*RP-DEBUG: pasted in from collisions.c:*/
			else {
				const double ls = 1.49597892e13;
				const double ts = 5.0226355648e6;
				double rn[3],vn=0,dTmp;
				for (i=0;i<3;i++) {
					rn[i] = c2.r[i] - c1.r[i];
					vn += (c2.v[i] - c1.v[i])*rn[i];
					}
				/* note vn >= 0 ==> near miss -- pkdBounce() will deal with it */
				if (vn < 0) {
					vn = - vn/sqrt(rn[0]*rn[0] + rn[1]*rn[1] + rn[2]*rn[2]);
					vn *= ls/ts; /* conversion to cm/s */
					switch (CP->iEpsNOption) {
					case PowerLaw:
						dEpsN = CP->dEpsNCoef*pow(vn,CP->dEpsNExp);
						break;
					case Compacted: /* Hatzes et al. 1988 */
						dTmp = c1.fRadius*c2.fRadius/(c1.fRadius + c2.fRadius)*ls;
						dEpsN = -1.90*exp(-(-0.01*dTmp + 0.41)*vn);
						break;
					case Borderies: /* Borderies et al. 1984 */
						if (vn >= CP->dEpsNVStar) {
							dTmp = CP->dEpsNVStar/vn;
							dTmp *= dTmp;
							dEpsN = sqrt(-2.0/3.0*dTmp + sqrt(10.0/3.0*dTmp - 5.0/9.0*dTmp*dTmp));
							}
						else
							dEpsN = 1.0;
						break;
					default:
						assert(0); /* Ran out of options */
						}
					if (dEpsN < CP->dEpsNMin)
						dEpsN = CP->dEpsNMin;
					else if (dEpsN > 1.0)
						dEpsN = 1.0;
					//dEpsT = CP->dEpsT; /* Not supported for AGGS..? */
					}
				}

#ifdef AGGS_IN_PATCH
			aggsBounce(&c1,&c2,dEpsN,&c,&n,pkd->PP,pkd->dTime + dt);
#else
			aggsBounce(&c1,&c2,dEpsN,&c,&n);
#endif			
			assert(n == 2);
#ifdef AGGS_IN_PATCH
			/* adjust y-momenta according to velocity and position (backdrift) change */
			/* Detailed explaination regarding Py [RP - 10/24/09]:
			  -At any time during the drift we can calculate the Py for any particle:
             Py = y_dot_n+1/2 + 2 Omega (x_n+? + x_dot_n+1/2 (dt/2 - t_n+?))
			   (cf. email with Tom Q, Sept 2009)

 			    Definitions:
			  -Py (cannonical momentum, similar to a specific angular momentum),
			  -y_dot, x_dot (velocities during the drift [time independent])
			  -Omega (orbital frequency)
			  -x_n+? (radial position of the particle at collision time)
			  -dt (timestep), t_n+? (elapsed time since the start of the drift)
			   ...the n+1/2 timeperiod is the 'drift.' n+? is the time of collision/event.
			   
			  -Any time a particle needs to know its Py, it uses this equation:  
			   It's used here, after a bounce; above, at fragmentation time; elsewhere,
			    when a particle breaks from an aggregate due to stress; and finally
			    when an aggregate is deleted due to only one particle remaining.
			  -It's used for COMs, as well as free particles.

			  -Aggregate COMs dictate the motion of the aggregate as a whole, so it is not
			    necessary for a constituent particle to keep track of its Py after every 
			    collision and timestep. However, when a particle is liberated, it is
			    important to update its Py at that moment, since its stored Py is
			    likely to be wrong at that later time.

			  -It is necessary to use this complex equation as Py is defined 
			    using the pre-kick velocity.  The velocities used during
			    the drift have already been kicked, and so this formula undoes that 
			    change.
			*/
			if (COLLIDER_IS_AGG(&c[0]) && c1.id.iPid == pkd->idSelf) {
				Aggregate *anew = &(c[0].agg);
				anew->dPy_com = anew->v_com[1] + 
				  2.0*pkd->PP->dOrbFreq*(anew->r_com[0] + anew->v_com[0]*(dHalfDelta - dt)); /*RP-DEBUG-dPy*/
			}

			if (COLLIDER_IS_AGG(&c[1]) && bReturnOutput)/* The first proc (c1's proc) handles all agg info */ {
				Aggregate *anew = &(c[1].agg);				
				anew->dPy_com = anew->v_com[1] + 
				  2.0*pkd->PP->dOrbFreq*(anew->r_com[0] + anew->v_com[0]*(dHalfDelta - dt)); /*RP-DEBUG-dPy*/
				/* RP-DEBUG: Remove BCs from ghosted particle as well as agg. */
				c[1].r[0] -= *dx;
				c[1].r[1] -= *dy;
				c[1].v[1] -= *dvy;
				c[1].dPy += *dvy/3.0;
				
				c[1].agg.r_com[0] -= *dx;
				c[1].agg.r_com[1] -= *dy;
				c[1].agg.v_com[1] -= *dvy;
				c[1].agg.dPy_com += *dvy/3.0;
			        }
#endif /*AGGS_IN_PATCH*/

			/*
			** Trace unaggregated particles back to start of step
			** (particles in aggregates updated in msrAggsBounce()).
			** The crazy (i+1)%2 business below is simply shorthand
			** for storing the iOrder of the *other* collider, since
			** for bouncing there can be only 2 particles involved.
			*/

			for (i=0;i<n;i++)
				if (!COLLIDER_IS_AGG(&c[i]) && c[i].id.iPid == pkd->idSelf) { 

#ifdef AGGS_IN_PATCH
				  /* RP:
				     Recalculate canonical y-momenta (Py) for unaggregated particles.
				     Yes, we can do this before backdrifting.  (Not a 
				     problem to do it afterward--we'd merely use another form of the 
				     equation.)
				  */
				        c[i].dPy = c[i].v[1] + 
					  2.0*pkd->PP->dOrbFreq*(c[i].r[0] + c[i].v[0]*(dHalfDelta - dt)); /*RP-DEBUG-dPy*/
#endif /* AGGS_IN_PATCH */					
					/* Backtrace */
				        for (k=0;k<3;k++)
						c[i].r[k] -= c[i].v[k]*dt;
#ifdef AGGS_IN_PATCH
					if(i==1) {
					        c[1].r[0] -= *dx;
						c[1].r[1] -= *dy;
						c[1].v[1] -= *dvy;
						c[1].dPy += *dvy/3.0;
					        }
#endif /* AGGS_IN_PATCH */
					pkdPutColliderInfo(pkd,&c[i],c[(i+1)%2].id.iOrder,&pkd->pStore[c[i].id.iIndex],dt);
				}
#ifdef AGGS_IN_PATCH
			*dx = *dy = *dvy = 0.0;
#endif

			/* cOut passes agg COM info back up to master. 
			   Only c1's aggregate info is reliable! */
			if (bReturnOutput) {
				for (i=0;i<n;i++) 
					cOut[i] = c[i]; /* struct copy */
				*pnOut = n;
				}
			/* free resources */
			free((void *) c);
			}
		}
	else {
		assert(0); /* no other outcomes allowed */
		}

	/*
	** For aggs, need to set dtPrevCol and reset iPrevCol.  (For
	** non-aggregate particles, this is done in pkdPutColliderInfo().)
	** Note we would like to set iPrevCol to INT_MAX (infinity) here
	** because it's possible for particles inside different aggs to
	** collide with one another multiple times during the interval
	** without intervening collisions involving those particles.
	** Unfortunately, possibly because of the approximate nature of
	** backdrifting, infinite loops have arisen due to repeated
	** collisions between the same aggregate (and agg/non-agg)
	** particles.  So, we must forbid these collisions, which means
	** it's possible to miss a legitimate collision.  Bottom line:
	** overlap correction is probably always necessary with
	** aggregates.  And small timesteps!
	*/

	if (COLLIDER_IS_AGG(&c1) && c1.id.iPid == pkd->idSelf) {
		pkd->pStore[c1.id.iIndex].dtPrevCol = dt;
		pkd->pStore[c1.id.iIndex].iPrevCol = c2.id.iOrder;
		}

	if (COLLIDER_IS_AGG(&c2) && c2.id.iPid == pkd->idSelf) {
		pkd->pStore[c2.id.iIndex].dtPrevCol = dt;
		pkd->pStore[c2.id.iIndex].iPrevCol = c1.id.iOrder;
		}
	}

/*** Following routines used for aggregate spin and orientation updates ***/

void aggsEulerDerivs(FLOAT t,FLOAT vars[],void *agg_as_void,FLOAT derivs[])
{
	/* t unused */

	FLOAT *torque = ((Aggregate *) agg_as_void)->torque;
	FLOAT *moments = ((Aggregate *) agg_as_void)->moments;

	/* omega[0, 1, 2] */
	derivs[0] = (torque[0] + vars[1]*vars[2]*(moments[1] - moments[2]))/moments[0];
	derivs[1] = (torque[1] + vars[2]*vars[0]*(moments[2] - moments[0]))/moments[1];
	derivs[2] = (torque[2] + vars[0]*vars[1]*(moments[0] - moments[1]))/moments[2];
 
	/* q1[0, 1, 2] */
	derivs[3] = vars[2]*vars[6] - vars[1]*vars[9];
	derivs[4] = vars[2]*vars[7] - vars[1]*vars[10];
	derivs[5] = vars[2]*vars[8] - vars[1]*vars[11];
 
	/* q2 */
	derivs[6] = vars[0]*vars[9]  - vars[2]*vars[3];
	derivs[7] = vars[0]*vars[10] - vars[2]*vars[4];
	derivs[8] = vars[0]*vars[11] - vars[2]*vars[5];
 
	/* q3 */
	derivs[9]  = vars[1]*vars[3] - vars[0]*vars[6];
	derivs[10] = vars[1]*vars[4] - vars[0]*vars[7];
	derivs[11] = vars[1]*vars[5] - vars[0]*vars[8];
	}

typedef void (*odeint_derivs)(FLOAT x,FLOAT y_vals[],void* user_data,
							  FLOAT deriv_vals[]);

/* following adapted from NRiC2e Ch. 16.2 (revised -- see online book) */

static
void rkck(FLOAT y[], FLOAT dydx[], int n, FLOAT x, FLOAT h, FLOAT yout[],
		  FLOAT yerr[], odeint_derivs derivs, void *user_data)
	 /*
	   Given values for n variables y[1..n] and their derivatives
	   dydx[1..n] known at x, use the fifth-order Cash-Karp
	   Runge-Kutta method to advance the solution over an interval h
	   and return the incremented variables as yout[1..n]. Also return
	   an estimate of the local truncation error in yout using the
	   embedded fourth-order method. The user supplies the routine
	   derivs(x,y,dydx), which returns derivatives dydx at x.
	   user_data is passed to derivs.
	 */
{
	int i;
	static FLOAT a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
	  b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
	  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
	  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
	  b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
	  c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
	  dc5 = -277.00/14336.0;
	FLOAT dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
	  dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	FLOAT *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
	ak2=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ak2 != NULL);
	ak3=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ak3 != NULL);
	ak4=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ak4 != NULL);
	ak5=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ak5 != NULL);
	ak6=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ak6 != NULL);
	ytemp=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ytemp != NULL);
	for (i=0;i<n;i++) /* First step. */
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,user_data,ak2); /* Second step. */
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,user_data,ak3); /* Third step. */
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,user_data,ak4); /* Fourth step. */
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,user_data,ak5); /* Fifth step. */
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,user_data,ak6); /* Sixth step. */
	for (i=0;i<n;i++) /* Accumulate increments with proper weights. */
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++) /* Estimate error as difference between fourth and fifth order methods. */
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	free((void *) ytemp);
	free((void *) ak6);
	free((void *) ak5);
	free((void *) ak4);
	free((void *) ak3);
	free((void *) ak2);
	}

#define SAFETY 0.9 /* SAFETY, PGROW, and PSHRNK relate to Eq. (16.2.10) */
#define PGROW  (-0.2)
#define PSHRNK (-0.25)
#define ERRCON 1.89e-4 /* ERRCON is (5/SAFETY)^(1/PGROW); see use below */

FLOAT FMAX(FLOAT,FLOAT);
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

FLOAT FMIN(FLOAT,FLOAT);
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static
void rkqs(FLOAT y[], FLOAT dydx[], int n, FLOAT *x, FLOAT htry, FLOAT eps,
		  FLOAT yscal[], FLOAT *hdid, FLOAT *hnext, odeint_derivs derivs,
		  void *user_data)
	 /*
	   Fifth-order Runge-Kutta step with monitoring of local
	   truncation error to ensure accuracy and adjust stepsize. Input
	   are the dependent variable vector y[1..n] and its derivative
	   dydx[1..n] at the starting value of the independent variable
	   x. Also input are the stepsize to be attempted htry, the
	   required accuracy eps, and the vector yscal[1..n] against which
	   the error is scaled. On output, y and x are replaced by their
	   new values, hdid is the stepsize that was actually
	   accomplished, and hnext is the estimated next stepsize. derivs
	   is the user-supplied routine that computes the right-hand side
	   derivatives. user_data is passed through to derivs.
	 */
{
	int i;
	FLOAT errmax,h,htemp,xnew,*yerr,*ytemp;
	FLOAT maxarg1,maxarg2,minarg1,minarg2; /* for FMAX() and FMIN() */
	yerr=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(yerr != NULL);
	ytemp=(FLOAT *) malloc(n*sizeof(FLOAT));
	assert(ytemp != NULL);
	h=htry; /* Set stepsize to the initial trial value. */
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs,user_data); /* Take a step. */
		errmax=0.0; /* Evaluate accuracy. */
		for (i=0;i<n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps; /* Scale relative to required tolerance. */
		if (errmax <= 1.0) break; /* Step succeeded. Compute size of next step. */
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		/* Truncation error too large, reduce stepsize. */
		h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h)); /* No more than a factor of 10. */
		xnew=(*x)+h;
		if (xnew == *x) {
			(void) fprintf(stderr,"rkqs(): stepsize underflow\n");
			assert(0);
			}
		}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h; /* No more than a factor of 5 increase. */
	*x += (*hdid=h);
	for (i=0;i<n;i++) y[i]=ytemp[i];
	free((void *) ytemp);
	free((void *) yerr);
	}

#undef FMIN
#undef FMAX
#undef ERRCON
#undef PSHRNK
#undef PGROW
#undef SAFETY

#define MAXSTP 10000
#define TINY 1.0e-30

FLOAT SIGN(FLOAT,FLOAT);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* code for storage of intermediate results omitted */

static
void odeint(FLOAT ystart[], int nvar, FLOAT x1, FLOAT x2, FLOAT eps, FLOAT h1,
			FLOAT hmin, int *nok, int *nbad, odeint_derivs derivs, void *user_data,
			void (*rkqs)(FLOAT [], FLOAT [], int, FLOAT *, FLOAT, FLOAT, FLOAT [],
						 FLOAT *, FLOAT *, odeint_derivs, void *))
	 /*
	   Runge-Kutta driver with adaptive stepsize control. Integrate
	   starting values ystart[1..nvar] from x1 to x2 with accuracy
	   eps. h1 should be set as a guessed first stepsize, hmin as the
	   minimum allowed stepsize (can be zero). On output nok and nbad
	   are the number of good and bad (but retried and fixed) steps
	   taken, and ystart is replaced by values at the end of the
	   integration interval. derivs is the user-supplied routine for
	   calculating the right-hand side derivative, while rkqs is the
	   name of the stepper routine to be used. user_data is passed
	   through to derivs.
	 */
{
	int nstp,i;
	FLOAT x,hnext,hdid,h;
	FLOAT *yscal,*y,*dydx;
	yscal=(FLOAT *) malloc(nvar*sizeof(FLOAT));
	assert(yscal != NULL);
	y=(FLOAT *) malloc(nvar*sizeof(FLOAT));
	assert(y != NULL);
	dydx=(FLOAT *) malloc(nvar*sizeof(FLOAT));
	assert(dydx != NULL);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = *nbad = 0;
	for (i=0;i<nvar;i++) y[i]=ystart[i];
	for (nstp=0;nstp<MAXSTP;nstp++) { /* Take at most MAXSTP steps. */
		(*derivs)(x,y,user_data,dydx);
		for (i=0;i<nvar;i++)
			/* Scaling used to monitor accuracy. This general-purpose
			   choice can be modified if need be. */
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x; /* If stepsize can overshoot, decrease. */
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs,user_data);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) { /* Are we done? */
			for (i=0;i<nvar;i++) ystart[i]=y[i];
			free((void *) dydx);
			free((void *) y);
			free((void *) yscal);
			return; /* Normal exit. */
			}
		if (fabs(hnext) <= hmin) {
			(void) fprintf(stderr,"WARNING: odeint(): small step size (%g)\n",fabs(hnext));
			/*assert(0);*/
			}
		h=hnext;
		assert(fabs(h) > 0.0);
		}
	(void) fprintf(stderr,"odeint(): too many steps\n");
	assert(0);
	}

#undef SIGN
#undef TINY
#undef MAXSTP

#define EPS 1.0e-6 /* desired fractional accuracy */
#define MIN_STEP 0.01 /* warn if step size < MIN_STEP * dDelta */
/*#define VERBOSE*/

void aggsRungeAdvance(Aggregate *agg,double dDelta,double dt)
{
	/* Interface to odeint() */

	FLOAT fVars[12],hmin;
	int k,nOK,nBad;

	/*
	** Build fVars.  To make things easier in aggsEulerDerivs(), we
	** copy the lambda data by *columns* (i.e. by principal axes).
	*/

	for (k=0;k<3;k++) {
		fVars[k] = agg->omega[k];
		fVars[k + 3] = agg->lambda[k][0]; /* q1 */
		fVars[k + 6] = agg->lambda[k][1]; /* q2 */
		fVars[k + 9] = agg->lambda[k][2]; /* q3 */
		}

	/* warn on small step, but only for long enough advance interval */

	hmin = (fabs(dt) <= MIN_STEP*dDelta ? 0.0 : MIN_STEP*dDelta);

	/* agg structure provides torque and moments */

	odeint(fVars,12,0.0,dt,EPS,dt,hmin,&nOK,&nBad,
		   aggsEulerDerivs,(void *) agg,rkqs);

#ifdef VERBOSE
	(void) printf("Runge advance from %g, interval %g: nOK = %i, nBad = %i\n",
				  agg->dLastUpdate,dt,nOK,nBad);
#endif

	/* Unpack fVars */

	for (k=0;k<3;k++) {
		agg->omega[k] = fVars[k];
		agg->lambda[k][0] = fVars[k + 3];
		agg->lambda[k][1] = fVars[k + 6];
		agg->lambda[k][2] = fVars[k + 9];
		}
	}

#undef VERBOSE
#undef EPS

/*#define TEST_RUNGE*/
#ifdef TEST_RUNGE
/*
** To test Runge-Kutta integrator, extract rkck(), rkqs(), odeint(),
** test_derivs(), and main() to separate file, add "#include
** <stdio.h>", "#include <math.h>", "#include <assert.h>", and
** "typedef double FLOAT;" at the top, compile (linking against the
** math library), and run.
*/

void test_derivs(FLOAT fDum,FLOAT fVars[],void *vDum,FLOAT derivs[])
{
	/* simple harmonic oscillator, amplitude 1, period 2 pi, phase 0 */

	derivs[0] = fVars[1];
	derivs[1] = - fVars[0];
	}

int main(void)
{
	FLOAT fVars[2],t,dt;
	int nOK,nBad;

	fVars[0] = 0.0;
	fVars[1] = 1.0;

	dt = 0.1;
	for (t=0.0;t<10.0;t+=dt) {
		odeint(fVars,2,t,t + dt,1.0e-9,dt,0.0,&nOK,&nBad,test_derivs,NULL,rkqs);
		(void) printf("t x v = %g %g %g (nOK = %i, nBad = %i)\n",t + dt,fVars[0],fVars[1],nOK,nBad);
		}

	return 0;
	}
#endif /* TEST_RUNGE */

#endif /* AGGS */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "mdl.h"
#include "pkd.h"
#include "grav.h"
#include "meval.h"
#include "qeval.h"

/*
 * Architectures for which the native sqrt() is faster.
 */

#define NATIVE_SQRT (defined(_MIPS_ISA) && (_MIPS_ISA == _MIPS_ISA_MIPS4) \
	     || defined(__i486__) || defined(CCC) || defined(__ia64__) \
	     || defined(__crayx1) || defined(OPTERON))

#ifdef NATIVE_SQRT
#undef NATIVE_SQRT
#endif

// This is always faster on modern hardware
#define NATIVE_SQRT 1

#if !(NATIVE_SQRT)
void v_sqrt1(int,double *,double *);
#endif

int pkdBucketInteract(PKD pkd,int iBucket,int iOrder)
{
    const KDN * const pkdn = &pkd->kdNodes[iBucket];
    PARTICLE * p = &pkd->pStore[pkdn->pLower];

    // get vals
    int nActive = 0;
    int n = pkdn->pUpper - pkdn->pLower + 1;
    int nPart = pkd->nPart;
    int nCellSoft = pkd->nCellSoft;
    int nCellNewt = pkd->nCellNewt;

#ifdef COMPLETE_LOCAL
    int nMultiFlop[5] = MEVAL_FLOP;
#else
    int nMultiFlop[5] = QEVAL_FLOP;
#endif

    const ILP * const ilp = pkd->ilp;
    const ILCS * const ilcs = pkd->ilcs;
    const ILCN * const ilcn = pkd->ilcn;
    int i;

	/*
	 ** Now process the two interaction lists for each active particle.
	 */
	for (i=0; i<n; ++i) {
        int j;
        double fPot,ax,ay,az;
        double x,y,z,dx,dy,dz,d2,h,twoh,a,b,c,d;
        double dir2,qirx,qiry,qirz,qir,tr,qir3;
        double idt2; /* reciprocal square of symmetric timestep */
        double gam[6];
        double dir;

		if (!TYPEQueryACTIVE(&(p[i])))
			continue;
		++nActive;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		fPot = 0.0;
		x = p[i].r[0];
		y = p[i].r[1];
		z = p[i].r[2];
		h = p[i].fSoft;

		/*
		 ** Scoring for Part (+,*)
		 ** 	Without sqrt = (10,8)
		 **     1/sqrt est.  = (6,11)
		 **     SPLINEM      = (0,3)  for soft = (8,30)
		 **     Total        = (16,22)           (24,49)
		 **     			 = 38	  for soft = 73
		 */
#if !(NATIVE_SQRT)
		for (j=0;j<nPart;++j) {
			dx = x - ilp[j].x;
			dy = y - ilp[j].y;
			dz = z - ilp[j].z;
			d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (nPart>0) v_sqrt1(nPart,d2a,sqrttmp);
#endif
		for (j=0; j<nPart; ++j) {
			dx = x - ilp[j].x;
			dy = y - ilp[j].y;
			dz = z - ilp[j].z;
			twoh = h + ilp[j].h;
#if (NATIVE_SQRT)
			d2 = dx*dx + dy*dy + dz*dz;
			SPLINE(d2,twoh,a,b);
#else
			SPLINEM(sqrttmp[j],d2a[j],twoh,a,b);
#endif
			idt2 = (p[i].fMass + ilp[j].m)*b;
			if (idt2 > p[i].dtGrav)
				p[i].dtGrav = idt2;
			a *= ilp[j].m;
			b *= ilp[j].m;
			fPot -= a;
			ax -= dx*b;
			ay -= dy*b;
			az -= dz*b;
			}

		/*
		 ** Scoring for CellSoft (+,*)
		 ** 	Without sqrt = (27,29)
		 **     1/sqrt est.  = (6,11)
		 **     SPLINEQ      = (0,9)  for soft = (13,62)
		 **     Total        = (33,49)           (46,102)
		 **     			 = 82	  for soft = 148
		 */
#if !(NATIVE_SQRT)
		for (j=0;j<nCellSoft;++j) {
			dx = x - ilcs[j].x;
			dy = y - ilcs[j].y;
			dz = z - ilcs[j].z;
			d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (nCellSoft>0) v_sqrt1(nCellSoft,d2a,sqrttmp);
#endif
		for (j=0;j<nCellSoft;++j) {
			dx = x - ilcs[j].x;
			dy = y - ilcs[j].y;
			dz = z - ilcs[j].z;
			twoh = h + ilcs[j].h;
#if (NATIVE_SQRT)
			d2 = dx*dx + dy*dy + dz*dz;
			dir = 1.0/sqrt(d2);
			SPLINEQ(dir,d2,twoh,a,b,c,d);
#else
			SPLINEQ(sqrttmp[j],d2a[j],twoh,a,b,c,d);
#endif
			qirx = ilcs[j].xx*dx + ilcs[j].xy*dy + ilcs[j].xz*dz;
			qiry = ilcs[j].xy*dx + ilcs[j].yy*dy + ilcs[j].yz*dz;
			qirz = ilcs[j].xz*dx + ilcs[j].yz*dy + ilcs[j].zz*dz;
			qir = 0.5*(qirx*dx + qiry*dy + qirz*dz);
			tr = 0.5*(ilcs[j].xx + ilcs[j].yy + ilcs[j].zz);
			qir3 = b*ilcs[j].m + d*qir - c*tr;
			fPot -= a*ilcs[j].m + c*qir - b*tr;
			ax -= qir3*dx - c*qirx;
			ay -= qir3*dy - c*qiry;
			az -= qir3*dz - c*qirz;
			idt2 = (p[i].fMass + ilcs[j].m)*b;
			if (idt2 > p[i].dtGrav)
				p[i].dtGrav = idt2;
			}
		/*
		 ** Try a cache check to improve responsiveness?
		 */
		mdlCacheCheck(pkd->mdl);
		/*
		 ** Scoring for CellNewt (+,*)
		 ** 	Without sqrt = (5,13)
		 **     1/sqrt est.  = (6,11)
		 **     Subtotal	 = (11,24)  = 35
		 **     Qeval        (Hex)		= 277 (see qeval.h)
		 **     Total        = (85,227) = 312 Flops/Newt-Interact
		 */
#if !(NATIVE_SQRT)
		for (j=0;j<nCellNewt;++j) {
			dx = x - ilcn[j].x;
			dy = y - ilcn[j].y;
			dz = z - ilcn[j].z;
			d2a[j] = dx*dx + dy*dy + dz*dz;
			}
		if (nCellNewt>0) v_sqrt1(nCellNewt,d2a,sqrttmp);
#endif
		for (j=0;j<nCellNewt;++j) {
			dx = x - ilcn[j].x;
			dy = y - ilcn[j].y;
			dz = z - ilcn[j].z;
#if (NATIVE_SQRT)
			d2 = dx*dx + dy*dy + dz*dz;
			gam[0] = 1.0/sqrt(d2);
#else
			gam[0] = sqrttmp[j];
#endif
			dir2 = gam[0]*gam[0];
			gam[1] = gam[0]*dir2;
			gam[2] = 3*gam[1]*dir2;
			gam[3] = 5*gam[2]*dir2;
			gam[4] = 7*gam[3]*dir2;
			gam[5] = 9*gam[4]*dir2;
#ifdef COMPLETE_LOCAL
			MEVAL(iOrder,ilcn[j],gam,dx,dy,dz,ax,ay,az,fPot);
#else
			QEVAL(iOrder,ilcn[j],gam,dx,dy,dz,ax,ay,az,fPot);  /* RP-DEBUG: acceleration alteration #3; source of voids! */
#endif
			idt2 = (p[i].fMass + ilcn[j].m)*gam[1];
			if (idt2 > p[i].dtGrav)
				p[i].dtGrav = idt2;
			}
		p[i].fPot += fPot;
		p[i].a[0] += ax;
		p[i].a[1] += ay;
		p[i].a[2] += az;
		/*
		 ** Try a cache check to improve responsiveness?
		 */
		mdlCacheCheck(pkd->mdl);
		}
	/*
	 ** Do the intra-bucket interactions.
	 ** Scoring (+,*):
	 ** 	without sqrt = (14,17)
	 **     sqrt est.    = (6,11)
	 **     SPLINE       = (0,3)  for soft = (8,30)
	 **     Total        = (20,31)           (28,58)
	 **                  = 51     for soft = 86
	 ** Multiplied by (n*(n-1)/2)!
	 */
	for (i=0;i<n-1;++i) {
		int j;

        double dx,dy,dz,twoh,a,b,d2;
        double idt2; /* reciprocal square of symmetric timestep */

#ifdef COLLISIONS
		double repelPrefactor; /* RP 6-22-09 */
#ifdef REPEL_MARK_II
		double dvx,dvy,dvz,dxDotDv;
		double totalMass;
#endif
#endif

		for (j=i+1;j<n;++j) {
			if (!TYPEQueryACTIVE(&(p[i])) && !TYPEQueryACTIVE(&(p[j])))
				continue;
			dx = p[j].r[0] - p[i].r[0];
			dy = p[j].r[1] - p[i].r[1];
			dz = p[j].r[2] - p[i].r[2];
			d2 = dx*dx + dy*dy + dz*dz;
			twoh = p[i].fSoft + p[j].fSoft;
#ifdef COLLISIONS
			/* RP 6-22-09: In case of repel-method overlap correction, compute linear repulsive force: */
			if (pkd->bRepel && d2 < (twoh)*(twoh)) {
				mdlassert(pkd->mdl,d2 > 0.0); /* Particles not allowed to perfectly overlap */

				/*
				** Code below, in REPEL_MARK_II, is commented pending
				** future improvements.  It is intended to remove
				** approach velocity from overlapping particles;
				** however, in practice, it only removes velocities
				** from UNAGGREGATED particles.  The edits to velocity
				** are ignored for particles in aggs, since an update
				** never occurs.  This is inconsistent, as a free
				** particle may overlap an agg, leading to a loss of
				** momentum conservation.  Perhaps this can be
				** resolved in the future.
				*/
#ifdef REPEL_MARK_II
				/*RP-DEBUG: (7/23/09) If particles are moving *toward* one another, nullify that motion */
				dvx = p[j].v[0] - p[i].v[0];
				dvy = p[j].v[1] - p[i].v[1];
				dvz = p[j].v[2] - p[i].v[2];
				
				dxDotDv = dx*dvx + dy*dvy + dz*dvz;
				if (dxDotDv < 0.0) /* That is, if particles are approaching one another */
				{
				    dxDotDv /= d2; /* For convenience */
				    totalMass = p[j].fMass+p[i].fMass;
				    
				    p[j].v[0] -= dx*dxDotDv*p[j].fMass/totalMass;
				    p[j].v[1] -= dy*dxDotDv*p[j].fMass/totalMass;
				    p[j].v[2] -= dz*dxDotDv*p[j].fMass/totalMass;

				    p[i].v[0] += dx*dxDotDv*p[i].fMass/totalMass;
				    p[i].v[1] += dy*dxDotDv*p[i].fMass/totalMass;
				    p[i].v[2] += dz*dxDotDv*p[i].fMass/totalMass;
				    
#ifdef SLIDING_PATCH
				    /* RP-DEBUG-dPy revision 11/5/09 */
				    /* Calculate Py based on repelled velocity */
				    /* NOTE: Uses a full timestep as the 'event time' */
				    p[i].dPy = p[i].v[1] + 2.0*pkd->PP->dOrbFreq*
				      (p[i].r[0] - (p[i].v[0]*pkd->PP->dDelta/2.0)); /*RP-DEBUG-dPy*/
				    p[j].dPy = p[j].v[1] + 2.0*pkd->PP->dOrbFreq*
				      (p[j].r[0] - (p[j].v[0]*pkd->PP->dDelta/2.0)); /*RP-DEBUG-dPy*/
#endif /* SLIDING_PATCH */
					}
#endif /* REPEL_MARK_II */
				
				/* Compute & apply repulsive force */
				/* Note: 1. pkd->dRepelFac is computed as 'k' in master.c: k = user_constant*dt^-2 */
				/*	     2. Using repulsive force law: F = -m*kx */		
				a = 1.0/sqrt(d2);
				b = -pkd->dRepelFac;
				repelPrefactor = (twoh*a) - 1.0;
				idt2 = -(p[i].fMass + p[j].fMass)*b; /* RP-DEBUG: b = -dt^-2 now (It has the same units as below.)*/
				if (TYPEQueryACTIVE(&(p[j]))) {
					p[j].fPot -= a*p[i].fMass;
					p[j].a[0] -= dx*repelPrefactor*b;
					p[j].a[1] -= dy*repelPrefactor*b;
					p[j].a[2] -= dz*repelPrefactor*b;
					if (idt2 > p[j].dtGrav)
						p[j].dtGrav = idt2;
					}
				if (TYPEQueryACTIVE(&(p[i]))) {
					p[i].fPot -= a*p[j].fMass;
					p[i].a[0] += dx*repelPrefactor*b;
					p[i].a[1] += dy*repelPrefactor*b;
					p[i].a[2] += dz*repelPrefactor*b;
					if (idt2 > p[i].dtGrav)
						p[i].dtGrav = idt2;
					}
				}
			else {
#endif /* COLLISIONS */
				SPLINE(d2,twoh,a,b); 
				idt2 = (p[i].fMass + p[j].fMass)*b;
				if (TYPEQueryACTIVE(&(p[j]))) {
					p[j].fPot -= a*p[i].fMass;
					p[j].a[0] -= dx*b*p[i].fMass;
					p[j].a[1] -= dy*b*p[i].fMass;
					p[j].a[2] -= dz*b*p[i].fMass;
					if (idt2 > p[j].dtGrav)
						p[j].dtGrav = idt2;
					}
				if (TYPEQueryACTIVE(&(p[i]))) {
					p[i].fPot -= a*p[j].fMass;
					p[i].a[0] += dx*b*p[j].fMass;
					p[i].a[1] += dy*b*p[j].fMass;
					p[i].a[2] += dz*b*p[j].fMass;
					if (idt2 > p[i].dtGrav)
						p[i].dtGrav = idt2;
					}
#ifdef COLLISIONS
				}
#endif
			}
		}
	/*
	 ** Compute the nFlop estimate.
	 */
	int nFlop = nActive*((nPart + n)*38 + nCellSoft*82 +
					 nCellNewt*(35 + nMultiFlop[iOrder]));
	return(nFlop);
	}

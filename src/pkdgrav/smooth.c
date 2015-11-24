#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "smoothfcn.h"


int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,
				 int bSymmetric,int iSmoothType,int bSmooth, double dfBall2OverSoft2 )
{
	SMX smx;
	void (*initParticle)(void *) = NULL;
	void (*initTreeParticle)(void *) = NULL;
	void (*init)(void *) = NULL;
	void (*comb)(void *,void *) = NULL;
	int pi;
	int nTree;

	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->pkd = smf->pkd = pkd;
	smx->nSmooth = nSmooth;
	smx->bPeriodic = bPeriodic;
	smx->dfBall2OverSoft2 = dfBall2OverSoft2;
	smx->iLowhFix = ( dfBall2OverSoft2 > 0.0 ? LOWHFIX_HOVERSOFT : LOWHFIX_NONE );
	smx->bUseBallMax = 1;
#ifdef SLIDING_PATCH
	/*DEBUG important: it is up to calling function to initialize these!*/
	smx->dTime = smf->dTime;
	smx->PP = &smf->PP;
#endif
#ifdef WALLS
	/*DEBUG important: it is up to calling function to initialize these!*/
	smx->WP = &smf->WP;
#endif

	switch (iSmoothType) {
	case SMX_NULL:
		smx->fcnSmooth = NullSmooth;
		initParticle = NULL; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = NULL; /* Cached copies */
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0;
		break;
	case SMX_DENSITY:
		smx->fcnSmooth = bSymmetric?DensitySym:Density;
		initParticle = initDensity; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = initDensity; /* Cached copies */
		comb = combDensity;
		smx->fcnPost = NULL;
		break;
	case SMX_MARKDENSITY:
		smx->fcnSmooth = bSymmetric?MarkDensitySym:MarkDensity;
		initParticle = initParticleMarkDensity; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = initMarkDensity; /* Cached copies */
		comb = combMarkDensity;
		smx->fcnPost = NULL;
		break;
	case SMX_MARKIIDENSITY:
		smx->fcnSmooth = bSymmetric?MarkIIDensitySym:MarkIIDensity;
		initParticle = initParticleMarkIIDensity; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = initMarkIIDensity; /* Cached copies */
		comb = combMarkIIDensity;
		smx->fcnPost = NULL;
		break;
	case SMX_MARK:
		smx->fcnSmooth = NULL;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = initMark;
		comb = combMark;
		smx->fcnPost = NULL;
		break;
	case SMX_DELTAACCEL:
		smx->fcnSmooth = DeltaAccel;
		initParticle = NULL; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = NULL; /* Cached copies */
		comb = combDeltaAccel;
		smx->fcnPost = NULL;
		break;
	case SMX_SINKACCRETE:
		smx->fcnSmooth = SinkAccrete;
		initParticle = NULL; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = initSinkAccrete; /* Cached copies */
		comb = combSinkAccrete;
		smx->fcnPost = NULL;
		smx->iLowhFix = LOWHFIX_SINKRADIUS;
		smx->bUseBallMax = 0;
		break;
	case SMX_BHSINKACCRETE:
		smx->fcnSmooth = BHSinkAccrete;
		initParticle = NULL; /* Original Particle */
		initTreeParticle = NULL; /* Original Particle */
		init = initBHSinkAccrete; /* Cached copies */
		comb = combBHSinkAccrete;
		smx->fcnPost = NULL;
		break;
#ifdef SUPERCOOL
	case SMX_MEANVEL:
		smx->fcnSmooth = bSymmetric?MeanVelSym:MeanVel;
		initParticle = initMeanVel;
		initTreeParticle = NULL;
		init = initMeanVel;
		comb = combMeanVel;
		smx->fcnPost = NULL;
		break;
#endif
#ifdef GASOLINE
	case SMX_SPHPRESSURETERMS:
		smx->fcnSmooth = bSymmetric?SphPressureTermsSym:SphPressureTerms;
		initParticle = initSphPressureTermsParticle; /* Original Particle */
		initTreeParticle = NULL;
		init = initSphPressureTerms; /* Cached copies */
		comb = combSphPressureTerms;
		smx->fcnPost = NULL;
		break;
	case SMX_DIVVORT:
		smx->fcnSmooth = bSymmetric?DivVortSym:DivVort;
		initParticle = initDivVort;
		initTreeParticle = NULL;
		init = initDivVort;
		comb = combDivVort;
		smx->fcnPost = NULL;
		break;
	case SMX_SHOCKTRACK:
		smx->fcnSmooth = bSymmetric?ShockTrackSym:ShockTrack;
		initParticle = initShockTrack;
		initTreeParticle = NULL;
		init = initShockTrack;
		comb = combShockTrack;
		smx->fcnPost = NULL;
		break;
	case SMX_HKPRESSURETERMS:
		smx->fcnSmooth = bSymmetric?HKPressureTermsSym:HKPressureTerms;
		initParticle = initHKPressureTermsParticle; /* Original Particle */
		initTreeParticle = NULL;
		init = initHKPressureTerms; /* Cached copies */
		comb = combHKPressureTerms;
		smx->fcnPost = NULL;
		break;
	case SMX_SPHPRESSURE:
		smx->fcnSmooth = bSymmetric?SphPressureSym:SphPressure;
		initParticle = initSphPressureParticle; /* Original Particle */
		init = initSphPressure; /* Cached copies */
		comb = combSphPressure;
		smx->fcnPost = postSphPressure;
		break;
	case SMX_SPHVISCOSITY:
		smx->fcnSmooth = bSymmetric?SphViscositySym:SphViscosity;
		initParticle = initSphViscosityParticle; /* Original Particle */
		initTreeParticle = NULL;
		init = initSphViscosity; /* Cached copies */
		comb = combSphViscosity;
		smx->fcnPost = NULL;
		break;
	case SMX_HKVISCOSITY:
		smx->fcnSmooth = bSymmetric?HKViscositySym:HKViscosity;
		initParticle = initHKViscosityParticle; /* Original Particle */
		initTreeParticle = NULL;
		init = initHKViscosity; /* Cached copies */
		comb = combHKViscosity;
		smx->fcnPost = NULL;
		break;
#ifdef STARFORM
	case SMX_DELETE_GAS:
                assert(bSymmetric == 0);
	        smx->fcnSmooth = DeleteGas;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0;
		break;
	case SMX_DIST_DELETED_GAS:
                assert(bSymmetric != 0);
	        smx->fcnSmooth = DistDeletedGas;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = initDistDeletedGas;
		comb = combDistDeletedGas;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0;
		break;
	case SMX_DIST_SN_ENERGY:
                assert(bSymmetric != 0);
	        smx->fcnSmooth = DistSNEnergy;
		initParticle = NULL;
		initTreeParticle = initTreeParticleDistSNEnergy;
		init = initDistSNEnergy;
		comb = combDistSNEnergy;
		smx->fcnPost = postDistSNEnergy;
		smx->bUseBallMax = 0;
		break;
#endif
#ifdef SIMPLESF
	case SMX_SIMPLESF_FEEDBACK:
		assert(bSymmetric != 0);
		smx->fcnSmooth = SimpleSF_Feedback;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = initSimpleSF_Feedback;
		comb = combSimpleSF_Feedback;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0;
		break;
#endif

#endif	       
#ifdef COLLISIONS
	case SMX_REJECTS:
		assert(bSymmetric != 0);
		smx->fcnSmooth = FindRejects;
		initParticle = initFindRejects;
		initTreeParticle = NULL;
		init = initFindRejects;
		comb = combFindRejects;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0;
		break;
	case SMX_COLLISION:
		assert(bSymmetric == 0);
		smx->fcnSmooth = CheckForCollision;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0; /* ignored */
		break;
#ifdef SPRINGS
	case SMX_ASSIGN_SPRINGS:
		assert(bSymmetric == 0); /* read-only cache */
		smx->fcnSmooth = AssignSprings;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0; /* ignored */
		break;
	case SMX_COLOR_SPRINGS:
		assert(bSymmetric == 0); /* read-only cache */
		smx->fcnSmooth = ColorSprings;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0; /* ignored */
		break;
	case SMX_MOVESPECIALS_SPRINGS:
		assert(bSymmetric == 0); /* read-only cache */
		smx->fcnSmooth = MoveSpecialSprings;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0; /* ignored */
		break;
	case SMX_DO_SPRINGS:
		assert(bSymmetric != 0);
		smx->fcnSmooth = DoSprings;
		initParticle = initDoSpringsParticle;
		initTreeParticle = NULL;
		init = initDoSprings;
		comb = combDoSprings;
		smx->fcnPost = postDoSprings;
		smx->bUseBallMax = 0; /* ignored */
		break;
#endif /* SPRINGS */
#ifdef DEM
	case SMX_ASSIGN_DEM:
		assert(bSymmetric == 0); /* read-only cache */
		smx->fcnSmooth = AssignDEM;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0; /* ignored */
		break;
	case SMX_DO_DEM:
		assert(bSymmetric != 0);
		smx->fcnSmooth = DoDEM;
		initParticle = initDoDEMParticle;
		initTreeParticle = NULL;
		init = initDoDEM;
		comb = combDoDEM;
		smx->fcnPost = postDoDEM;
		smx->bUseBallMax = 0; /* ignored */
		break;
	case SMX_DEM_STATS:
		assert(bSymmetric == 0); /* read-only cache */
		smx->fcnSmooth = DEMStats;
		initParticle = NULL;
		initTreeParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		smx->bUseBallMax = 0; /* ignored */
		break;
#endif /* DEM */
#ifdef RORY_EXTENSION
	case SMX_FINDBINARY:
		smx->fcnSmooth = FindBinary;
		initParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		break;
#endif /* RORY_EXTENSION */
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
	case SMX_FIND_OVERLAPS:
		assert(bSymmetric != 0);
		smx->fcnSmooth = FindOverlaps;
		initParticle = initFindOverlaps;
		initTreeParticle = NULL;
		init = initFindOverlaps;
		comb = combFindOverlaps;
		smx->fcnPost = NULL;
		break;
#endif /* SLIDING_PATCH */
	default:
		assert(0);
		}
	/*
	 ** Initialize the ACTIVE particles in the tree.
	 ** There are other particles in the tree -- just not active.
	 */
	nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
	for (pi=0;pi<nTree;++pi) {
		if (TYPEQuerySMOOTHACTIVE(&(pkd->pStore[pi]))) {
			TYPEReset( &(pkd->pStore[pi]),TYPE_SMOOTHDONE );
			/*			if (bSmooth) pkd->pStore[pi].fBall2 = -1.0;*/
			if (initParticle != NULL) {
				initParticle(&pkd->pStore[pi]);
				}
			}
                else if (TYPEQueryTREEACTIVE(&(pkd->pStore[pi]))) {
			/*TYPEReset( &(pkd->pStore[pi]),TYPE_SMOOTHDONE );*/
			/*			if (bSmooth) pkd->pStore[pi].fBall2 = -1.0;*/
			if (initTreeParticle != NULL) {
				initTreeParticle(&pkd->pStore[pi]);
				}
			}
		}
	/*
	 ** Start particle caching space (cell cache is already active).
	 */
	if (bSymmetric) {
		assert(init != NULL);
		assert(comb != NULL);
		mdlCOcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
				   nTree,init,comb);
		}
	else {
		mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
				   nTree);
		}
	/*
	 ** Allocate mark array.
	 */
	if(pkdLocal(pkd) == 0)
	    smx->piMark = (int *)malloc(sizeof(int));
	else
	    smx->piMark = (int *)malloc(pkdLocal(pkd)*sizeof(int));
	assert(smx->piMark != NULL);
	/*
	 ** Allocate Nearest-Neighbor List.
	 */
	smx->nListSize = smx->nSmooth;
	smx->nnList = (NN *)malloc(smx->nListSize*sizeof(NN));
	assert(smx->nnList != NULL);
	smx->pbRelease = (int *)malloc(smx->nListSize*sizeof(int));
	assert(smx->pbRelease != NULL);
	/*
	 ** Allocate priority queue.
	 */
	smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
	assert(smx->pq != NULL);
	PQ_INIT(smx->pq,nSmooth);
	/*
	 ** Set up "in priority queue" hash table stuff.
	 */
	smx->nHash = (int)ceil(nSmooth/PQ_LOAD_FACTOR);
	smx->pqHash = (PQ **)malloc(smx->nHash*sizeof(PQ *));
	assert(smx->pqHash != NULL);
	*psmx = smx;	
	return(1);
	}


void smFinish(SMX smx,SMF *smf, CASTAT *pcs)
{
	PKD pkd = smx->pkd;
	int pi;

	/*
	 ** Get caching statistics.
	 */
	pcs->dcNumAccess = mdlNumAccess(pkd->mdl,CID_CELL);
	pcs->dcMissRatio = mdlMissRatio(pkd->mdl,CID_CELL);
	pcs->dcCollRatio = mdlCollRatio(pkd->mdl,CID_CELL);
	pcs->dcMinRatio = mdlMinRatio(pkd->mdl,CID_CELL);
	pcs->dpNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE);
	pcs->dpMissRatio = mdlMissRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpCollRatio = mdlCollRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpMinRatio = mdlMinRatio(pkd->mdl,CID_PARTICLE);
	/*
	 ** Stop particle caching space.
	 */
	mdlFinishCache(smx->pkd->mdl,CID_PARTICLE);
	/*
	 ** Now do any post calculations, these ususlly involve some sort of
	 ** normalizations of the smoothed quantities, usually division by
	 ** the local density! Do NOT put kernel normalizations in here as
	 ** these do not depend purely on local properties in the case of
	 ** "Gather-Scatter" kernel.
	 */
	if (smx->fcnPost != NULL) {
		for (pi=0;pi<pkd->nTreeActive;++pi) {
                        smx->fcnPost(&pkd->pStore[pi],smf);	
			}
		}
	/*
	 ** Free up context storage.
	 */
	free(smx->pqHash);
	free(smx->pq);
	free(smx->piMark);
	free(smx->nnList);
	free(smx->pbRelease);
	free(smx);
	}


PQ *smBallSearch(SMX smx,PQ *pq,FLOAT *ri,int *cpStart)
{
	MDL mdl = smx->pkd->mdl;
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int cell,cp,ct,pj,pUpper,idSelf,bPeriodic;
	FLOAT fBall2,fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	if(smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLOAT_MAXVAL;
	    ly = FLOAT_MAXVAL;
	    lz = FLOAT_MAXVAL;
	    }
	idSelf = smx->pkd->idSelf;
	fBall2 = pq->fKey;
	cell = ROOT;
	bPeriodic = 0;
	/*
	 ** First find the "local" Bucket.
	 ** This could be the closest to ri[3]! 
	 ** Warning: in parallel ri[3] SHOULD be contained in the LOCAL DOMAIN!
	 */
	while (c[cell].iDim >= 0) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	pUpper = c[cell].pUpper;
	for (pj=c[cell].pLower;pj<=pUpper;++pj) {
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			*cpStart = cell;
			if (smx->piMark[pj]) continue;
			if (pq->id == idSelf) smx->piMark[pq->p] = 0;
			else {
				mdlRelease(mdl,CID_PARTICLE,pq->pPart);
				PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
				pq->id = idSelf;
				}
			smx->piMark[pj] = 1;
			pq->fKey = fDist2;
			pq->dx = dx;
			pq->dy = dy;
			pq->dz = dz;
			pq->p = pj;
			pq->pPart = &p[pj];
			pq->ax = 0.0;
			pq->ay = 0.0;
			pq->az = 0.0;
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,
					  GetNext_1);
			if (c[cp].iDim >= 0) {
				cp = LOWER(cp);
				continue;
				}
			else {
				pUpper = c[cp].pUpper;
				for (pj=c[cp].pLower;pj<=pUpper;++pj) {
					dx = sx - p[pj].r[0];
					dy = sy - p[pj].r[1];
					dz = sz - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						*cpStart = PARENT(cell);
						if (smx->piMark[pj]) continue;
						if (pq->id == idSelf) smx->piMark[pq->p] = 0;
						else {
							mdlRelease(mdl,CID_PARTICLE,pq->pPart);
							PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
							pq->id = idSelf;
							}
						smx->piMark[pj] = 1;
						pq->fKey = fDist2;
						pq->dx = dx;
						pq->dy = dy;
						pq->dz = dz;
						pq->p = pj;
						pq->pPart = &p[pj];
						pq->ax = sx - x;
						pq->ay = sy - y;
						pq->az = sz - z;
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNext_1:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	if (bPeriodic) *cpStart = 0;
	return(pq);
	}


PQ *smBallSearchNP(SMX smx,PQ *pq,FLOAT *ri,int *cpStart)
{
	MDL mdl = smx->pkd->mdl;
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int cell,cp,ct,pj,pUpper,idSelf;
	FLOAT fBall2,fDist2,dx,dy,dz,x,y,z;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	idSelf = smx->pkd->idSelf;
	fBall2 = pq->fKey;
	cell = ROOT;
	/*
	 ** First find the "local" Bucket.
	 ** This could mearly be the closest to ri[3]! 
	 ** Warning: in parallel ri[3] SHOULD be contained in the LOCAL DOMAIN!
	 */
	while (c[cell].iDim >= 0) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	pUpper = c[cell].pUpper;
	for (pj=c[cell].pLower;pj<=pUpper;++pj) {
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			*cpStart = cell;

			if (smx->piMark[pj])
                continue;

			if (pq->id == idSelf)
            {
                smx->piMark[pq->p] = 0;
            }
			else
            {
				mdlRelease(mdl,CID_PARTICLE,pq->pPart);
				PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
				pq->id = idSelf;
            }

			smx->piMark[pj] = 1;
			pq->fKey = fDist2;
			pq->dx = dx;
			pq->dy = dy;
			pq->dz = dz;
			pq->p = pj;
			pq->pPart = &p[pj];
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECTNP(&c[cp],fBall2,x,y,z,GetNext_1);
			if (c[cp].iDim >= 0) {
				cp = LOWER(cp);
				continue;
				}
			else {
				pUpper = c[cp].pUpper;
				for (pj=c[cp].pLower;pj<=pUpper;++pj) {
					dx = x - p[pj].r[0];
					dy = y - p[pj].r[1];
					dz = z - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						*cpStart = PARENT(cell);
						if (smx->piMark[pj]) continue;
						if (pq->id == idSelf) smx->piMark[pq->p] = 0;
						else {
							mdlRelease(mdl,CID_PARTICLE,pq->pPart);
							PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
							pq->id = idSelf;
                        }
						smx->piMark[pj] = 1;
						pq->fKey = fDist2;
						pq->dx = dx;
						pq->dy = dy;
						pq->dz = dz;
						pq->p = pj;
						pq->pPart = &p[pj];
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNext_1:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	return(pq);
	}


/*
 * I don't increase the PQ size here: I assume it isn't in use if we
 * need this function.
 */
void smGrowList(SMX smx)
{
    smx->nListSize *= 1.5;
    
    smx->nnList = (NN *) realloc(smx->nnList,smx->nListSize*sizeof(NN));
    assert(smx->nnList != NULL);
    smx->pbRelease =
		(int *) realloc(smx->pbRelease,smx->nListSize*sizeof(int));
    assert(smx->pbRelease != NULL);
}

int smBallGather(SMX smx,FLOAT fBall2,FLOAT *ri)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int pj,nCnt,cp,pUpper;
	FLOAT dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
	int iDum;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	if (smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLOAT_MAXVAL;
	    ly = FLOAT_MAXVAL;
	    lz = FLOAT_MAXVAL;
	    }
	nCnt = 0;
	cp = ROOT;
	if(c[cp].pLower > c[cp].pUpper) {
	    return(0);		/* Empty Tree */
	    }
	    
	while (1) {
		INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = sx - p[pj].r[0];
				dy = sy - p[pj].r[1];
				dz = sz - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= fBall2) {
					if(nCnt >= smx->nListSize)
					    smGrowList(smx);
					smx->nnList[nCnt].fDist2 = fDist2;
					smx->nnList[nCnt].dx = dx;
					smx->nnList[nCnt].dy = dy;
					smx->nnList[nCnt].dz = dz;
					smx->nnList[nCnt].pPart = &p[pj];
					smx->nnList[nCnt].iIndex = pj;
					smx->nnList[nCnt].iPid = smx->pkd->idSelf;
					smx->pbRelease[nCnt++] = 0;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}


int smBallGatherNP(SMX smx,FLOAT fBall2,FLOAT *ri,int cp)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int pj,nCnt,pUpper;
	FLOAT dx,dy,dz,x,y,z,fDist2;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	nCnt = 0;
	if(c[cp].pLower > c[cp].pUpper) {
	    return(0);		/* Empty Tree */
	    }
	while (1) {
		INTERSECTNP(&c[cp],fBall2,x,y,z,GetNextCell);
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = x - p[pj].r[0];
				dy = y - p[pj].r[1];
				dz = z - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= fBall2) {
					if(nCnt >= smx->nListSize)
					    smGrowList(smx);
					smx->nnList[nCnt].fDist2 = fDist2;
					smx->nnList[nCnt].dx = dx;
					smx->nnList[nCnt].dy = dy;
					smx->nnList[nCnt].dz = dz;
					smx->nnList[nCnt].pPart = &p[pj];
					smx->nnList[nCnt].iIndex = pj;
					smx->nnList[nCnt].iPid = smx->pkd->idSelf;
					smx->pbRelease[nCnt++] = 0;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}


void smBallScatter(SMX smx,FLOAT *ri,int iMarkType)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int pj,cp,pUpper;
	FLOAT dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
	int iDum;

	/*mdlDiag(smx->pkd->mdl, "smBallScatter:Start\n" );*/
	x = ri[0];
	y = ri[1];
	z = ri[2];
	if (smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLOAT_MAXVAL;
	    ly = FLOAT_MAXVAL;
	    lz = FLOAT_MAXVAL;
	    }

	cp = ROOT;
	while (1) {
		/*mdlDiag(smx->pkd->mdl, "smINTERSECTSCATTER:Before\n" );*/
		INTERSECTSCATTER(&c[cp],lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
		/*mdlDiag(smx->pkd->mdl, "smINTERSECTSCATTER:After\n" );*/
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = sx - p[pj].r[0];
				dy = sy - p[pj].r[1];
				dz = sz - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= p[pj].fBallMax*p[pj].fBallMax) {
					TYPESet(&(p[pj]),iMarkType);
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}

	/*mdlDiag(smx->pkd->mdl, "smBallScatter:End\n" );*/
	}


void smBallScatterNP(SMX smx,FLOAT *ri,int iMarkType,int cp)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int pj,pUpper;
	FLOAT dx,dy,dz,x,y,z,fDist2;

	x = ri[0];
	y = ri[1];
	z = ri[2];

	while (1) {
		INTERSECTSCATTERNP(&c[cp],x,y,z,GetNextCell);
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = x - p[pj].r[0];
				dy = y - p[pj].r[1];
				dz = z - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= p[pj].fBallMax*p[pj].fBallMax) {
					TYPESet(&(p[pj]),iMarkType);
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	}


void smSmooth(SMX smx,SMF *smf)
{
	PKD pkd = smx->pkd;
	MDL mdl = smx->pkd->mdl;
	KDN *c = pkd->kdNodes;
	PARTICLE *p = pkd->pStore;
	KDN *pkdn;
	PARTICLE *pPart;
	int nSmooth,i,j,pi,pj,pNext,nCnt;
	int cell,idcell,cp,id,ct,idct;
	FLOAT fBall2,fDist2,x,y,z,dx,dy,dz,lx,ly,lz,sx,sy,sz,h2;
	PQ *pq,*pqi,*pqn;
	int iDum;
	int nTree,nQueue,iLoad;
	int nLocal;
	int nSmoothed = 0;
	int idSelf;
	KDN dpkdn; /* dummy node */

	idSelf = smx->pkd->idSelf;
	nSmooth = smx->nSmooth;
	nTree = c[pkd->iRoot].pUpper + 1;
	nLocal = pkd->nLocal;
	if (smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLOAT_MAXVAL;
	    ly = FLOAT_MAXVAL;
	    lz = FLOAT_MAXVAL;
	    }
	/*
	 ** Clear Mark array and pqHash.
	 */
	for (pi=0;pi<pkd->nLocal;++pi) smx->piMark[pi] = 0;
	for (j=0;j<smx->nHash;++j) smx->pqHash[j] = NULL;
	pNext = 0;
 StartParticle:
	/*
	 ** Check if we are finished!
	 */
	if (pNext == nLocal) {
		goto DoneSmooth;
		}

	if (!TYPEFilter(&(pkd->pStore[pNext]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,TYPE_SMOOTHACTIVE)) {
		/*) || pkd->pStore[pNext].fBall2 >= 0.0) {*/
		++pNext;
		goto StartParticle;
		}
	pi = pNext;
	x = p[pi].r[0];
	y = p[pi].r[1];
	z = p[pi].r[2];
	cell = ROOT;
	while (c[cell].iDim >= 0) {
		if (p[pi].r[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	/*
	 ** Add local stuff to the prioq.
	 */
	iLoad = c[cell].pLower;
	nQueue = 0;
	if (iLoad > nTree - nSmooth) iLoad = nTree - nSmooth;
	/*
	 ** Have to add non-local particles to the queue?
	 */
	for (id=0;id<mdlThreads(pkd->mdl) && (iLoad < 0);++id) {
		if (id == pkd->idSelf) continue;
		pkdn = mdlAquire(mdl,CID_CELL,ROOT,id);
		for (pj=pkdn->pLower;pj<=pkdn->pUpper && (iLoad < 0);++pj) {
			pqi = &smx->pq[nQueue];
			pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
			dx = x - pPart->r[0];
			dy = y - pPart->r[1];
			dz = z - pPart->r[2];
			pqi->ax = 0.0;
			pqi->ay = 0.0;
			pqi->az = 0.0;
			/*
			 ** If nSmooth ~ N, may need to load the queue with ghosts...
			 */
			for (i=0;i<3;i++)
				dpkdn.bnd.fMin[i] = dpkdn.bnd.fMax[i] = pPart->r[i];
			INTERSECT(&dpkdn,FLOAT_MAXVAL,lx,ly,lz,x,y,z,sx,sy,sz,iDum,
					  dumlabel1);
			dx = sx - pPart->r[0];
			dy = sy - pPart->r[1];
			dz = sz - pPart->r[2];
			pqi->ax = sx - x;
			pqi->ay = sy - y;
			pqi->az = sz - z;
		dumlabel1:
			pqi->fKey = dx*dx + dy*dy + dz*dz;
			pqi->dx = dx;
			pqi->dy = dy;
			pqi->dz = dz;
			pqi->p = pj;
			pqi->id = id;
			pqi->pPart = pPart;
			PQ_HASHADD(smx->pqHash,smx->nHash,pqi);
			/*
			 ** Note: in this case we DO NOT want to release pPart!
			 */
			++iLoad;
			++nQueue;
			}
		mdlRelease(mdl,CID_CELL,pkdn);
		}
	assert(iLoad >= 0);
	for (;nQueue<nSmooth;++nQueue,++iLoad) {
		pqi = &smx->pq[nQueue];
		smx->piMark[iLoad] = 1;
		dx = x - p[iLoad].r[0];
		dy = y - p[iLoad].r[1];
		dz = z - p[iLoad].r[2];
		pqi->ax = 0.0;
		pqi->ay = 0.0;
		pqi->az = 0.0;
		for (i=0;i<3;i++)
			dpkdn.bnd.fMin[i] = dpkdn.bnd.fMax[i] = p[iLoad].r[i];
		INTERSECT(&dpkdn,FLOAT_MAXVAL,lx,ly,lz,x,y,z,sx,sy,sz,iDum,dumlabel2);
		dx = sx - p[iLoad].r[0];
		dy = sy - p[iLoad].r[1];
		dz = sz - p[iLoad].r[2];
		pqi->ax = sx - x;
		pqi->ay = sy - y;
		pqi->az = sz - z;
	dumlabel2:
		pqi->fKey = dx*dx + dy*dy + dz*dz;
		pqi->dx = dx;
		pqi->dy = dy;
		pqi->dz = dz;
		pqi->p = iLoad;
		pqi->id = pkd->idSelf;
		pqi->pPart = &p[iLoad];
		}
	PQ_BUILD(smx->pq,nSmooth,pq);
/*
	sprintf(ach,"Start:%d\n",p[pi].iOrder);
	mdlDiag(pkd->mdl,ach);
*/
 PrioqBuilt:
	/*
	 ** Priority Queue must be built. 'pi' must be defined.
	 */
	if (smx->bPeriodic) {
		pq = smBallSearch(smx,pq,p[pi].r,&p[pi].cpStart);
		}
	else {
		pq = smBallSearchNP(smx,pq,p[pi].r,&p[pi].cpStart);
	    }
	/*
	 ** Start non-local search.
	 */
    fBall2 = pq->fKey;
	idcell = -1;	/* We are in the LTT now ! */
	cell = pkd->piLeaf[pkd->idSelf];
	while (!pkdIsRoot(cell,idcell)) {
		cp = SIBLING(cell);
		id = idcell;
		ct = cp;
		idct = id;
		pkdNext(pkd,ct,idct);
		/*
		 ** Check for any cache work.
		 */
		while (1) {
		    if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
			else pkdn = &pkd->kdTop[cp];
			if (pkdn->pUpper < 0) goto GetNext_2;
   			INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNext_2);
			if (pkdn->iDim >= 0 || id == -1) {
				if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
				pkdLower(pkd,cp,id);
				continue;
				}
			for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
				pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
				dx = sx - pPart->r[0];
				dy = sy - pPart->r[1];
				dz = sz - pPart->r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 < fBall2) {
					p[pi].cpStart = 0;
					PQ_INQUEUE(smx->pqHash,smx->nHash,pj,id,NextParticle);
					if (pq->id == pkd->idSelf) smx->piMark[pq->p] = 0;
					else {
						mdlRelease(mdl,CID_PARTICLE,pq->pPart);
						PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
						}
					pq->fKey = fDist2;
					pq->dx = dx;
					pq->dy = dy;
					pq->dz = dz;
					pq->p = pj;
					pq->id = id;
					pq->pPart = pPart;
					pq->ax = sx - x;
					pq->ay = sy - y;
					pq->az = sz - z;
					PQ_HASHADD(smx->pqHash,smx->nHash,pq);
					PQ_REPLACE(pq);
					fBall2 = pq->fKey;
					/*
					 ** Note: in this case we DO NOT want to release pPart!
					 */
					/*
					 ** Setting p[pi].cpStart to zero means we have to search 
					 ** non-local tree when doing a BallGather (smReSmooth).
					 */
					continue;
					}
			NextParticle:
				mdlRelease(mdl,CID_PARTICLE,pPart);
				}
		GetNext_2:
			if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
			pkdNext(pkd,cp,id);
			if (cp == ct && id == idct) break;
			}
		pkdParent(pkd,cell,idcell);
		}
	/*
	 ** Create Nearest-Neighbor List and try to pick next particle.
	 */
	pqn = NULL;
	nCnt = 0;
	h2 = 2.0*pq->fKey;   /* arbitrarily bigger than pq->fKey! */

	fBall2 = pq->fKey;
	/* 
	** If desired reject fBall2 below a minimum 
	** We may get many more than nSmooth neighbours here -- resort to a ReSmooth
	*/
	if (smx->iLowhFix != LOWHFIX_NONE && 
		((smx->iLowhFix==LOWHFIX_HOVERSOFT && fBall2 < smx->dfBall2OverSoft2*p[pi].fSoft*p[pi].fSoft) ||
 		 (smx->iLowhFix==LOWHFIX_SINKRADIUS && fBall2 < smf->dSinkRadius*smf->dSinkRadius)) ) {
		/* We ReSmooth for this guy later */
		p[pi].fBall2 = -1.0; /* any value < 0 will do -- see code after "DoneSmooth:" below */
		TYPESet(&p[pi],TYPE_SMOOTHDONE);

		/* Get the next in line */
		for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
			if (pqi->id != pkd->idSelf) continue;
			if (!TYPEFilter(&(pkd->pStore[pqi->p]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,TYPE_SMOOTHACTIVE)) continue;
			if (pqi->fKey < h2) {
				pqn = pqi;
				h2 = pqn->fKey;
				}
			}
		}
	else {
		/* Limit fBall2 growth to help stability and neighbour finding */
		if (smx->bUseBallMax && p[pi].fBallMax > 0.0 && fBall2 > p[pi].fBallMax*p[pi].fBallMax)
			fBall2=p[pi].fBallMax*p[pi].fBallMax;

		p[pi].fBall2 = fBall2;
		TYPESet(&p[pi],TYPE_SMOOTHDONE);

		for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
			/*
			 ** There are cases like in collisions where we do want to include the
			 ** single nearest neighbor. So include the most distant particle.
			 **
			 if (pqi == pq) continue;
			 */
			/*
			 ** Move relevant data into Nearest Neighbor array.
			 */
			if (pqi->fKey <= fBall2) {
				smx->nnList[nCnt].iPid = pqi->id;
				smx->nnList[nCnt].iIndex = pqi->p;
				smx->nnList[nCnt].pPart = pqi->pPart;
				smx->nnList[nCnt].fDist2 = pqi->fKey;
				smx->nnList[nCnt].dx = pqi->dx;
				smx->nnList[nCnt].dy = pqi->dy;
				smx->nnList[nCnt].dz = pqi->dz;
				++nCnt;
				} 
			if (pqi->id != pkd->idSelf) continue;
			if (!TYPEFilter(&(pkd->pStore[pqi->p]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,TYPE_SMOOTHACTIVE)) continue;
			if (pqi->fKey < h2) {
				pqn = pqi;
				h2 = pqn->fKey;
				}
			}

#if (INTERNAL_WARNINGS != 0)

		/*
		** For periodic boundary conditions, it's probably a good idea
		** that the search ball not exceed half the spatial period.
		** However, there are enough special circumstances that this
		** is now always only a warning, not an error.
		*/

		if (smx->bPeriodic &&
			((lx < FLOAT_MAXVAL && p[pi].fBall2 >= 0.25*lx*lx) ||
			 (ly < FLOAT_MAXVAL && p[pi].fBall2 >= 0.25*ly*ly) ||
			 (lz < FLOAT_MAXVAL && p[pi].fBall2 >= 0.25*lz*lz))) {
			static int nWarn = 1;
			if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
				fprintf(stderr,"LARGE SEARCH BALL WARNING #%i (pid=%i): %i (fBall = %g; x,y,z = %g,%g,%g)\n",
						nWarn,pkd->idSelf,p[pi].iOrder,sqrt(p[pi].fBall2),p[pi].r[0],p[pi].r[1],p[pi].r[2]);
			++nWarn;
		}

#endif

		nSmoothed++;
		smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
	}

	/*
	 ** Need to do a CACHE recombine (ie. finish up Smooth)
	 ** to deliver info to particles on other processors.
	 */
	/*
	 ** Try a cache check to improve responsiveness.
	 */
	mdlCacheCheck(mdl);
	if (!pqn) {
		/*
		 ** Clean up piMark array, pqHash, and release all aquired
		 ** pointers!
		 */
		for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
			if (pqi->id == pkd->idSelf) smx->piMark[pqi->p] = 0;
			else {
				mdlRelease(mdl,CID_PARTICLE,pqi->pPart);
				PQ_HASHDEL(smx->pqHash,smx->nHash,pqi);
				}
			}
		goto StartParticle;
		}
	/*
	 ** Calculate the priority queue using the previous particles.
	 ** "THE SNAKE"
	 */
	pi = pqn->p;
	x = p[pi].r[0];
	y = p[pi].r[1];
	z = p[pi].r[2];
	for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
		if (pqi == pqn) continue;
		pqi->ax -= pqn->ax;
		pqi->ay -= pqn->ay;
		pqi->az -= pqn->az;
		dx = x - pqi->pPart->r[0] + pqi->ax;
		dy = y - pqi->pPart->r[1] + pqi->ay;
		dz = z - pqi->pPart->r[2] + pqi->az;
		pqi->fKey = dx*dx + dy*dy + dz*dz;
		pqi->dx = dx;
		pqi->dy = dy;
		pqi->dz = dz;
		}
	pqn->fKey = 0.0;
	pqn->ax = 0.0;
	pqn->ay = 0.0;
	pqn->az = 0.0;
	PQ_BUILD(smx->pq,nSmooth,pq);
/*
	sprintf(ach,"Snake:%d\n",p[pi].iOrder);
	mdlDiag(pkd->mdl,ach);
*/
	goto PrioqBuilt;

 DoneSmooth:
	for (pi=0;pi<nLocal;++pi) {
	    if (!TYPETest(&(p[pi]),TYPE_SMOOTHACTIVE) || p[pi].fBall2 >= 0.0) continue;
		/*
		 ** Do a Ball Gather at r^2 = fBall2
		 */
		switch (smx->iLowhFix) {
		case LOWHFIX_HOVERSOFT:
		  fBall2 = smx->dfBall2OverSoft2*p[pi].fSoft*p[pi].fSoft;
		  break;
		case LOWHFIX_SINKRADIUS:
		  fBall2 = smf->dSinkRadius*smf->dSinkRadius;
		  break;
		default:
		  fprintf(stderr,"Illegal value for iLowhFix %d in smooth (active = %i, fBall2 = %g)\n",smx->iLowhFix,TYPETest(&(p[pi]),TYPE_SMOOTHACTIVE),p[pi].fBall2);
		  assert(0);
		}
		p[pi].fBall2 = fBall2;
		p[pi].cpStart = 0;
		cp = 0;
		if (cp) {
			nCnt = smBallGatherNP(smx,fBall2,p[pi].r,cp);
			smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
			/*
			** Try a cache check to improve responsiveness.
			*/
			mdlCacheCheck(mdl);
			}
		else {
			if (smx->bPeriodic) {
				nCnt = smBallGather(smx,fBall2,p[pi].r);
				}
			else {
				nCnt = smBallGatherNP(smx,fBall2,p[pi].r,ROOT);
				}
			/*
			 ** Start non-local search.
			 */
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			cp = ROOT;
			id = -1;	/* We are in the LTT now ! */
			while (1) {
				if (id == pkd->idSelf) goto SkipLocal;
				if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
				else pkdn = &pkd->kdTop[cp];
				if (pkdn->pUpper < 0) goto GetNextCell;
				INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
				if (pkdn->iDim >= 0 || id == -1) {
					if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
					pkdLower(pkd,cp,id);
					continue;
					}
				for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
					pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
					dx = sx - pPart->r[0];
					dy = sy - pPart->r[1];
					dz = sz - pPart->r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 <= fBall2) {
						if(nCnt >= smx->nListSize)
						    smGrowList(smx);
						smx->nnList[nCnt].fDist2 = fDist2;
						smx->nnList[nCnt].dx = dx;
						smx->nnList[nCnt].dy = dy;
						smx->nnList[nCnt].dz = dz;
						smx->nnList[nCnt].pPart = pPart;
						smx->nnList[nCnt].iIndex = pj;
						smx->nnList[nCnt].iPid = id;
						smx->pbRelease[nCnt++] = 1;
						continue;
						}
					mdlRelease(mdl,CID_PARTICLE,pPart);
					}
			GetNextCell:
				if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
			SkipLocal:
				pkdNext(pkd,cp,id);
				if (pkdIsRoot(cp,id)) break;
				}

			smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
			/*
			 ** Release non-local particle pointers.
			 */
			for (i=0;i<nCnt;++i) {
				if (smx->pbRelease[i]) {
					mdlRelease(mdl,CID_PARTICLE,smx->nnList[i].pPart);
					}
				}
			}
		}
	/*
	printf("nSmoothed: %d\n",nSmoothed);
	*/
	}

void smReSmooth(SMX smx,SMF *smf)
{
	PKD pkd = smx->pkd;
	MDL mdl = smx->pkd->mdl;
	PARTICLE *p = smx->pkd->pStore;
	PARTICLE *pPart;
	KDN *pkdn;
	int pi,pj,nCnt,cp,id,i;
	FLOAT x,y,z,lx,ly,lz,sx,sy,sz,dx,dy,dz,fDist2,fBall2;
	int iDum;
	int nTree;
	
	if (smx->bPeriodic) {
	    lx = pkd->fPeriod[0];
	    ly = pkd->fPeriod[1];
	    lz = pkd->fPeriod[2];
	    }
	else {
	    lx = FLOAT_MAXVAL;
	    ly = FLOAT_MAXVAL;
	    lz = FLOAT_MAXVAL;
	    }
	nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
	for (pi=0;pi<nTree;++pi) {
		if (!TYPEFilter(&(p[pi]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,
						TYPE_SMOOTHACTIVE)) continue;

		/*
		 ** Do a Ball Gather at the radius of the most distant particle
		 ** which smSmooth sets in p[pi].fBall2.
		 */
		fBall2 = p[pi].fBall2;
		cp = p[pi].cpStart;
		if (nTree == 1 && fBall2 == 0.0)
			cp = 0; /* DCR: under these circumstances, cp may have nonsensical value */
		TYPESet(&p[pi],TYPE_SMOOTHDONE);
		if (cp) {
			nCnt = smBallGatherNP(smx,fBall2,p[pi].r,cp);
			smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
			/*
			** Try a cache check to improve responsiveness.
			*/
			mdlCacheCheck(mdl);
			}
		else {
			if (smx->bPeriodic) {
				nCnt = smBallGather(smx,fBall2,p[pi].r);
				}
			else {
				nCnt = smBallGatherNP(smx,fBall2,p[pi].r,ROOT);
				}
			/*
			 ** Start non-local search.
			 */
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			cp = ROOT;
			id = -1;	/* We are in the LTT now ! */
			while (1) {
				if (id == pkd->idSelf) goto SkipLocal;
				if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
				else pkdn = &pkd->kdTop[cp];
				if (pkdn->pUpper < 0) goto GetNextCell;
				INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
				if (pkdn->iDim >= 0 || id == -1) {
					if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
					pkdLower(pkd,cp,id);
					continue;
					}
				for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
					pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
					dx = sx - pPart->r[0];
					dy = sy - pPart->r[1];
					dz = sz - pPart->r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 <= fBall2) {
						if(nCnt >= smx->nListSize)
						    smGrowList(smx);
						smx->nnList[nCnt].fDist2 = fDist2;
						smx->nnList[nCnt].dx = dx;
						smx->nnList[nCnt].dy = dy;
						smx->nnList[nCnt].dz = dz;
						smx->nnList[nCnt].pPart = pPart;
						smx->nnList[nCnt].iIndex = pj;
						smx->nnList[nCnt].iPid = id;
						smx->pbRelease[nCnt++] = 1;
						continue;
						}
					mdlRelease(mdl,CID_PARTICLE,pPart);
					}
			GetNextCell:
				if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
			SkipLocal:
				pkdNext(pkd,cp,id);
				if (pkdIsRoot(cp,id)) break;
				}
			smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
			/*
			 ** Release non-local particle pointers.
			 */
			for (i=0;i<nCnt;++i) {
				if (smx->pbRelease[i]) {
					mdlRelease(mdl,CID_PARTICLE,smx->nnList[i].pPart);
					}
				}
			}
		}
	}

void smMarkSmooth(SMX smx,SMF *smf,int iMarkType)
{
	PKD pkd = smx->pkd;
	MDL mdl = smx->pkd->mdl;
	PARTICLE *p = smx->pkd->pStore;
	PARTICLE *pPart;
	KDN *pkdn;
	int pi,pj,cp,id;
	FLOAT x,y,z,lx,ly,lz,sx,sy,sz,dx,dy,dz,fDist2;
	int iDum;
	int nTree;
	char DiagStr[100];
	
	mdlDiag(smx->pkd->mdl,"smMarkSmooth:Start\n" );
	if (smx->bPeriodic) {
	    lx = pkd->fPeriod[0];
	    ly = pkd->fPeriod[1];
	    lz = pkd->fPeriod[2];
	    }
	else {
	    lx = FLOAT_MAXVAL;
	    ly = FLOAT_MAXVAL;
	    lz = FLOAT_MAXVAL;
	    }
	nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
	for (pi=0;pi<nTree;++pi) {
		if (!TYPEFilter(&(p[pi]),TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE,
						TYPE_SMOOTHACTIVE)) continue;
		/*		sprintf( DiagStr, "smMarkSmooth: Particle %d\n", p[pi].iOrder);*/
 		mdlDiag(smx->pkd->mdl, DiagStr );
		/*
		 ** Do a Ball Scatter to this particle
		 ** which is smSmooth sets in p[pi].fBall2.
		 */
		TYPESet(&p[pi],TYPE_SMOOTHDONE);

		if (smx->bPeriodic) {
			smBallScatter(smx,p[pi].r,iMarkType);
			}
		else {
			smBallScatterNP(smx,p[pi].r,iMarkType,ROOT);
			}
		/*
		 ** Start non-local search.
		 */
		x = p[pi].r[0];
		y = p[pi].r[1];
		z = p[pi].r[2];
		cp = ROOT;
		id = -1;	/* We are in the LTT now ! */
		while (1) {
			if (id == pkd->idSelf) goto SkipLocal;
			if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
			else pkdn = &pkd->kdTop[cp];
			if (pkdn->pUpper < 0) goto GetNextCell;
			INTERSECTSCATTER(pkdn,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
			if (pkdn->iDim >= 0 || id == -1) {
				if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
				pkdLower(pkd,cp,id);
				continue;
				}
			for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
				pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
				dx = sx - pPart->r[0];
				dy = sy - pPart->r[1];
				dz = sz - pPart->r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= pPart->fBallMax*pPart->fBallMax) 
				  TYPESet(pPart, iMarkType);	
				mdlRelease(mdl,CID_PARTICLE,pPart);
				}
		GetNextCell:
			if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
		SkipLocal:
			pkdNext(pkd,cp,id);
			if (pkdIsRoot(cp,id)) break;
			}
		}
	}

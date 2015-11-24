/*
 ** rpu.c -- DCR 98-09-24
 ** =====
 ** Generic utilities for rubble pile data.
 */

#include <rpu.h>
#include <stdlib.h>
#ifdef sgi
#include <bstring.h>
#else
#include <string.h>
#endif
#ifdef sparc
#include <strings.h>
#endif
#include <assert.h>
#include <boolean.h>

static void jacobi(MATRIX a00, VECTOR d0, MATRIX v00);
static void invert(MATRIX a);

double rpuVolSph(double r)
{
	return 4.0/3*PI*CUBE(r);
	}

double rpuVolEll(const VECTOR a)
{
	return 4.0/3*PI*a[X]*a[Y]*a[Z];
	}

BOOLEAN rpuInEllipsoid(const VECTOR r0,double R,const VECTOR a)
{
	/*
	 ** Returns TRUE if a ball with radius R located at r0 lies *entirely*
	 ** within the ellipsoid defined by semi-axes "a" (measured along the
	 ** Cartesian axes and centred at the origin). To get this right, need
	 ** to compute direction from r0 to nearest point on ellipsoid surface.
	 ** This is too hard, so settle for more conservative boundary.
	 */

	VECTOR r;
	double d;
	int k;

	assert(a[X] > 0.0 && a[Y] > 0.0 && a[Z] > 0.0);

	/* Check easy cases first */

	d = 0;

	for (k=0;k<N_DIM;k++) {
		if (r0[k] == 0 && R > a[k]) return FALSE;
		d += SQ(r0[k]/a[k]);
		}

	if (d == 0) return TRUE;
	if (d > 1) return FALSE;

	/* Make a (conservative) guess for the general case */

	d = MAG(r0);
	COPY_VEC(r0,r);
	SCALE_VEC(r,1 + R/d);
	d = SQ(r[X]/a[X]) + SQ(r[Y]/a[Y]) + SQ(r[Z]/a[Z]);

	return (d <= 1);
	}

/*
 ** Following functions take rubble-pile arguments.
 */

void rpuMalloc(RUBBLE_PILE *rp)
{
	assert(rp != NULL);
	assert(rp->n_particles > 0);
	rp->data = (SSDATA *) malloc(rp->n_particles*sizeof(SSDATA));
	assert(rp->data != NULL);
	}

void rpuRealloc(RUBBLE_PILE *rp)
{
	/* note: ok if rp->data is NULL ... equivalent to rpuMalloc(rp) */

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	rp->data =
		(SSDATA *) realloc((void *) rp->data,rp->n_particles*sizeof(SSDATA));
	assert(rp->data != NULL);
	}

void rpuFree(RUBBLE_PILE *rp)
{
	assert(rp != NULL);
	assert(rp->data != NULL);
	free((void *) rp->data);
	rp->data = NULL;
	}

void rpuCopyData(const RUBBLE_PILE *rps,RUBBLE_PILE *rpd,const size_t offset)
{
	bcopy((void *) rps->data,(void *) &rpd->data[offset],
		  rps->n_particles*sizeof(SSDATA));
	}

void rpuCalcMass(RUBBLE_PILE *rp)
{
	/* Calculates total mass of rubble pile */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	rp->mass = 0;
	for (i=0;i<rp->n_particles;i++)
		rp->mass += rp->data[i].mass;
	}

void rpuCalcRadius(RUBBLE_PILE *rp)
{
	/*
	 ** Estimates bulk radius of rubble pile (needs center of mass
	 ** position -- see rpuCalcPos()).
	 */

	VECTOR r;
	double r_mag;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	rp->radius = 0;
	for (i=0;i<rp->n_particles;i++) {
		SUB_VEC(rp->data[i].pos,rp->pos,r);
		r_mag = MAG(r) + rp->data[i].radius;
		if (r_mag > rp->radius) rp->radius = r_mag;
		}
	}

void rpuCalcPos(RUBBLE_PILE *rp)
{
	/*
	 ** Sets rp->pos to center-of-mass position. Note rp->mass must be
	 ** the total mass. Call rpuCalcMass() to compute this if necessary.
	 */

	VECTOR r;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	assert(rp->mass >= 0.0); /* zero mass now allowed */

	ZERO_VEC(rp->pos);
	for (i=0;i<rp->n_particles;i++) {
		COPY_VEC(rp->data[i].pos,r);
		SCALE_VEC(r,rp->data[i].mass);
		ADD_VEC(rp->pos,r,rp->pos);
		}
	if (rp->mass > 0.0) {
		NORM_VEC(rp->pos,rp->mass);
		}
	else {
		NORM_VEC(rp->pos,rp->n_particles);
		}
	}

void rpuCalcVel(RUBBLE_PILE *rp)
{
	/* Sets rp->vel to center-of-mass velocity -- see rpuCalcPos() */

	VECTOR v;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	assert(rp->mass >= 0.0); /* zero mass now allowed */

	ZERO_VEC(rp->vel);
	for (i=0;i<rp->n_particles;i++) {
		COPY_VEC(rp->data[i].vel,v);
		SCALE_VEC(v,rp->data[i].mass);
		ADD_VEC(rp->vel,v,rp->vel);
		}
	if (rp->mass > 0.0) {
		NORM_VEC(rp->vel,rp->mass);
		}
	else {
		NORM_VEC(rp->vel,rp->n_particles);
		}
	}

void rpuCalcInertia(RUBBLE_PILE *rp)
{
	/* Calculates inertia tensor of rubble pile (needs rp->pos) */

	/*
	** NOTE: The inertia tensor consists of the point mass
	** distribution PLUS the component spheres.  This allows us then
	** to set L = Iw to solve for the bulk spin w (on the assumption
	** that each component sphere has the same spin w).
	*/

	VECTOR r;
	int i,j,k;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	ZERO_MAT(rp->inertia);
	for (i=0;i<rp->n_particles;i++) {
		SUB_VEC(rp->data[i].pos,rp->pos,r);
		for (j=0;j<N_DIM;j++)
			for (k=0;k<N_DIM;k++)
				rp->inertia[j][k] +=
					rp->data[i].mass*((j == k ? 0.4*SQ(rp->data[i].radius) +
									   MAG_SQ(r) : 0) - r[j]*r[k]);
		}
	}

void rpuCalcAxes(RUBBLE_PILE *rp)
{
	/*
	 ** Calculates moments (eigenvalues) and axes (eigenvectors) of
	 ** rp->inertia. The principal axes are the ROWS of rp->axes.
	 ** This routine automatically calls rpuCalcAxisLengths().
	 */

	MATRIX m;
	double x,y,z;
	int k;

	assert(rp != NULL);

	COPY_MAT(rp->inertia,m); /* jacobi() modifies the matrix */
	jacobi(m,rp->moments,rp->axes);
	Transpose(rp->axes);
	for (k=0;k<N_DIM;k++) { /* make unit vectors */
		x = MAG(rp->axes[k]);
		assert(x > 0.0);
		NORM_VEC(rp->axes[k],x);
		}

	/* Determine moment order (smallest to largest) */

	x = rp->moments[X]; y = rp->moments[Y]; z = rp->moments[Z];

	if (x < y && x < z)
		if (y < z) SET_VEC(rp->mom_ord,X,Y,Z)
		else SET_VEC(rp->mom_ord,X,Z,Y)
	else if (y < z)
		if (x < z) SET_VEC(rp->mom_ord,Y,X,Z)
		else SET_VEC(rp->mom_ord,Y,Z,X)
	else
		if (x < y) SET_VEC(rp->mom_ord,Z,X,Y)
		else SET_VEC(rp->mom_ord,Z,Y,X)

	rpuCalcAxisLengths(rp);
	}

void rpuCalcAxisLengths(RUBBLE_PILE *rp)
{
	/*
	 ** Calculates the lengths of the rubble pile principal semi-axes. The
	 ** matrix rp->axes must contain the (normalized) axis orientations
	 ** (see rpuCalcAxes()). Optionally, rp->axes can be set to the unit
	 ** matrix to get body dimensions along the Cartesian axes.
	 ** NOTE: It is very hard to calculate the size of the best-fitting
	 ** ellipsoid exactly, so we use an iterative procedure to converge
	 ** on a reasonable guess. Note the lengths include the finite size of
	 ** the particles and assume spatial symmetry wrt the center of mass.
	 */

	/*DEBUG! updated to use center of figure, not center of mass, and to NOT
	  use rpuInEllipsoid()*/

	SSDATA *d;
	VECTOR r,s,vMin,vMax,cof;
	double x,y,z;
	int i,k;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	/* get maximum dimensions along body axes */

	ZERO_VEC(rp->axis_len);
	SET_VEC(rp->axis_ord,X,Y,Z);
	ZERO_VEC(vMin);
	ZERO_VEC(vMax);
	for (i=0;i<rp->n_particles;i++) {
		d = &rp->data[i];
		SUB_VEC(d->pos,rp->pos,r);
		Transform(rp->axes,r,s);
		for (k=0;k<N_DIM;k++) {
/*DEBUG! old...
			a = fabs(s[k]) + d->radius;
			if (a > rp->axis_len[k]) rp->axis_len[k] = a;*/
			if (s[k] - d->radius < vMin[k]) vMin[k] = s[k] - d->radius;
			if (s[k] + d->radius > vMax[k]) vMax[k] = s[k] + d->radius;
			}
		}

	/*DEBUG! new version*/

	for (k=0;k<N_DIM;k++) {
		rp->axis_len[k] = 0.5*(vMax[k] - vMin[k]); /* semi-axes */
		cof[k] = 0.5*(vMin[k] + vMax[k]);/*DEBUG not used!*/
		}
	goto order;

	/* expand uniformly in small increments until all particles contained within ellipsoid */

	for (i=0;i<rp->n_particles;i++) {
		d = &rp->data[i];
		SUB_VEC(d->pos,rp->pos,r);
		Transform(rp->axes,r,s);
		while (!rpuInEllipsoid(s,d->radius,rp->axis_len)) {
			SCALE_VEC(rp->axis_len,1.001); /* lengths good to 0.1% as N --> inf */
			}
		}

	/* determine axis order (longest to shortest) */

 order:/*DEBUG!*/
	x = rp->axis_len[X]; y = rp->axis_len[Y]; z = rp->axis_len[Z];

	if (x > y && x > z)
		if (y > z) SET_VEC(rp->axis_ord,X,Y,Z)
		else SET_VEC(rp->axis_ord,X,Z,Y)
	else if (y > z)
		if (x > z) SET_VEC(rp->axis_ord,Y,X,Z)
		else SET_VEC(rp->axis_ord,Y,Z,X)
	else
		if (x > y) SET_VEC(rp->axis_ord,Z,X,Y)
		else SET_VEC(rp->axis_ord,Z,Y,X)
	}

#define TOLERANCE (1.0e-12)

void rpuCalcSpin(RUBBLE_PILE *rp)
{
	/*
	 ** Calculates spin (and angular momentum) of rubble pile. Note the
	 ** inertia tensor and axis moments are required in advance -- see
	 ** rpuCalcInertia() and rpuCalcAxes().
	 */

	SSDATA *d;
	MATRIX a;
	VECTOR r,v,h;
	double Isph,Ix,Iy,Iz,Id;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	assert(rp->mass >= 0.0);

	ZERO_VEC(rp->ang_mom);
	if (rp->mass == 0.0) {
		ZERO_VEC(rp->spin);
		}
	else if (rp->n_particles == 1) {
		d = rp->data;
		Isph = 0.4*d->mass*SQ(d->radius);
		COPY_VEC(d->spin,rp->ang_mom);
		SCALE_VEC(rp->ang_mom,Isph);
		UNIT_MAT(rp->inertia);
		SCALE_MAT(rp->inertia,Isph);
		COPY_VEC(d->spin,rp->spin);
		}
	else {
		for (i=0;i<rp->n_particles;i++) {
			d = &rp->data[i];
			Isph = 0.4*d->mass*SQ(d->radius);
			SUB_VEC(d->pos,rp->pos,r);
			SUB_VEC(d->vel,rp->vel,v);
			CROSS(r,v,h);
			SCALE_VEC(h,d->mass);
			ADD_VEC(rp->ang_mom,h,rp->ang_mom);
			COPY_VEC(d->spin,h);
			SCALE_VEC(h,Isph);
			ADD_VEC(rp->ang_mom,h,rp->ang_mom);
			}
		COPY_MAT(rp->inertia,a);
		invert(a);
		Transform(a,rp->ang_mom,rp->spin);
		}

	/* Now compute rotation state parameters (Cf. Scheeres et al. 2000) */

	rp->ang_mom_mag = MAG(rp->ang_mom);
	rp->kin_energy = 0.5*DOT(rp->spin,rp->ang_mom);
	if (rp->ang_mom_mag == 0.0)
		rp->eff_spin = 0.0;
	else
		rp->eff_spin = 2*rp->kin_energy/rp->ang_mom_mag;
	if (rp->kin_energy == 0.0)
		rp->dyn_inertia = 0.0;
	else
		rp->dyn_inertia = 0.5*SQ(rp->ang_mom_mag)/rp->kin_energy;
	Ix = rp->moments[rp->mom_ord[X]]; /* min */
	Iy = rp->moments[rp->mom_ord[Y]]; /* mid */
	Iz = rp->moments[rp->mom_ord[Z]]; /* max */
	assert(Ix <= Iy && Iy <= Iz);
	Id = rp->dyn_inertia;
	if (Id == 0.0) {
		rp->rot_idx = 0.0;
		return;
		}
	if (Id < Ix && fabs(Id - Ix) <= TOLERANCE) Id = Ix;
	if (Id > Iz && fabs(Id - Iz) <= TOLERANCE) Id = Iz;
	assert(Id >= Ix && Id <= Iz);
	if (Ix == Iy && Ix == Iz) { /* spherical case */
		rp->rot_idx = 1.0;
		}
	else if (Ix == Iy && Ix < Iz) { /* oblate case */
		rp->rot_idx = (Id - Ix)/(Iz - Ix);
		}
	else if (Iz == Iy && Ix < Iz) { /* prolate case */
		rp->rot_idx = (Id - Iz)/(Iz - Ix);
		}
	else if (Id >= Iy) {
		rp->rot_idx = (Id - Iy)/(Iz - Iy);
		}
	else {
		rp->rot_idx = (Id - Iy)/(Iy - Ix);
		}
	if (rp->rot_idx > 1.0 && fabs(rp->rot_idx - 1.0) <= TOLERANCE) rp->rot_idx = 1.0;
	if (fabs(rp->rot_idx) <= TOLERANCE) rp->rot_idx = 0.0;
	if (rp->rot_idx < -1.0 && fabs(rp->rot_idx + 1.0) <= TOLERANCE) rp->rot_idx = -1.0;
	assert(rp->rot_idx >= -1.0 && rp->rot_idx <= 1.0);
	}

#undef TOLERANCE

void rpuCalcColor(RUBBLE_PILE *rp)
{
	/* Sets rubble pile color to most dominant particle color */

	int n[NUM_COLORS],c,n_max,i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	for (i=0;i<NUM_COLORS;i++)
		n[i] = 0;

	for (i=0;i<rp->n_particles;i++) {
		c = rp->data[i].color;
		if (c < 0)
			c = 0;
		else if (c >= NUM_COLORS)
			c = NUM_COLORS - 1;
		++n[c];
		}

	rp->color = n_max = 0;
	for (i=0;i<NUM_COLORS;i++)
		if (n[i] > n_max) {
			rp->color = i;
			n_max = n[i];
			}

	assert(rp->color >= 0 && rp->color < NUM_COLORS);
	}

void rpuCalcAggID(RUBBLE_PILE *rp)
{
	/*
	 ** Determines whether all particles in rubble pile share same
	 ** (negative) original index (org_idx).  If so, the aggregate
	 ** ID is computed as -1 - org_idx; otherwise it is set to -1.
	 */

	int i,agg_id;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	rp->agg_id = -1;
	for (i=0;i<rp->n_particles;i++) {
		agg_id = AGG_IDX(&rp->data[i]);
		if (agg_id < 0)
			return;
		if (agg_id != rp->agg_id) {
			if (rp->agg_id < 0)
				rp->agg_id = agg_id;
			else {
				rp->agg_id = -1;
				return;
				}
			}
		}
	}

void rpuCalcDensity(RUBBLE_PILE *rp)
{
	/*
	 ** Calculates (ellipsoidal) bulk density of rubble pile. The axis
	 ** lengths and mass must be known in advance.
	 */

	assert(rp != NULL);
	rp->density = rp->mass/rpuVolEll(rp->axis_len);
	}

void rpuCalcPacking(RUBBLE_PILE *rp)
{
	/*
	 ** Calculates particle packing efficiency. The particle data and
	 ** rubble pile axis lengths are required in advance.
	 */

	double volume = 0;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	for (i=0;i<rp->n_particles;i++)
		volume += rpuVolSph(rp->data[i].radius);

	rp->packing = volume/rpuVolEll(rp->axis_len);
	}

void rpuScaleMass(RUBBLE_PILE *rp,double f)
{
	/*
	 ** Scales particle mass and total mass by "f" while keeping particle
	 ** radii and bulk dimension constant. It is up to the calling routine
	 ** to recompute any derived quantities (e.g. inertia tensor, etc.).
	 */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	assert(f >= 0.0); /* zero mass now allowed */

	for (i=0;i<rp->n_particles;i++)
		rp->data[i].mass *= f;

	rp->mass *= f;
	}

void rpuScaleRadius(RUBBLE_PILE *rp,double f,BOOLEAN just_particles)
{
	/*
	 ** Scales particle radius and, if !just_particles, bulk
	 ** dimensions by "f" while keeping particle masses and total mass
	 ** constant.  It is up to the calling routine to recompute any
	 ** derived quantities (e.g. axis lengths).  Note that particle
	 ** positions are scaled out as well (if !just_particles).  Also,
	 ** if only the particles radii are being changed, the calling
	 ** program should recompute the bulk mass...
	 */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	if (just_particles && f < 0.0) {
		/* -ve value means set every particle radius to -f */
		for (i=0;i<rp->n_particles;i++)
			rp->data[i].radius = -f;
		return;
		}

	assert(f > 0.0);

	for (i=0;i<rp->n_particles;i++)
		rp->data[i].radius *= f;

	if (just_particles)
		return;

	for (i=0;i<rp->n_particles;i++)
		SCALE_VEC(rp->data[i].pos,f);

	rp->radius *= f;
	}

void rpuApplyPos(RUBBLE_PILE *rp)
{
	/* Sets new particle positions relative to rp->pos */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	for (i=0;i<rp->n_particles;i++)
		ADD_VEC(rp->data[i].pos,rp->pos,rp->data[i].pos);
	}

void rpuApplyVel(RUBBLE_PILE *rp)
{
	/* Sets new particle velocities relative to rp->vel */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	for (i=0;i<rp->n_particles;i++)
		ADD_VEC(rp->data[i].vel,rp->vel,rp->data[i].vel);
	}

void rpuApplyColor(RUBBLE_PILE *rp)
{
	/* Applies rubble pile color to particles */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	for (i=0;i<rp->n_particles;i++)
		rp->data[i].color = rp->color;
	}

void rpuApplyAggID(RUBBLE_PILE *rp)
{
	/* Assigns particles to single aggregate */

	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	for (i=0;i<rp->n_particles;i++)
		if (rp->agg_id >= 0) {
			AGG_SET_IDX(&rp->data[i],rp->agg_id);
			}
		else
			rp->data[i].org_idx = i;
	}

void rpuRotate(RUBBLE_PILE *rp,MATRIX rot0,BOOLEAN body)
{
	/*
	 ** Rotates rubble pile. If "body" is TRUE, use body axes, else space
	 ** axes.  After the transformation, rpuAnalyze() is called to
	 ** recompute other affected quantities, e.g. inertia tensor, etc.
	 */

	SSDATA *d;
	MATRIX rot,m;
	VECTOR v;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	COPY_MAT(rot0,rot); /* to ensure passed matrix is unchanged */

	if (body) {
		COPY_MAT(rot,m);
		MultMat(rp->axes,m,rot);
		COPY_MAT(rot,m);
		Transpose(rp->axes);
		MultMat(m,rp->axes,rot);
		Transpose(rp->axes);
		}
	for (i=0;i<rp->n_particles;i++) {
		d = &rp->data[i];
		SUB_VEC(d->pos,rp->pos,v);
		Transform(rot,v,d->pos);
		ADD_VEC(d->pos,rp->pos,d->pos);
		SUB_VEC(d->vel,rp->vel,v);
		Transform(rot,v,d->vel);
		ADD_VEC(d->vel,rp->vel,d->vel);
		COPY_VEC(d->spin,v);
		Transform(rot,v,d->spin);
		}
	rpuAnalyze(rp);
	}

void rpuAddSpin(RUBBLE_PILE *rp,BOOLEAN body)
{
	/*
	 ** Adds spin to particle velocities, taking rubble pile
	 ** position and velocity into account.  If "body" is TRUE,
	 ** added spin is applied relative to body axes, else space
	 ** axes.  To completely reset the spin, call this function
	 ** with rp->spin set to (- rp->spin) first---this will set
	 ** the bulk spin to zero.
	 */

	SSDATA *d;
	MATRIX m;
	VECTOR v,r;
	int i;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	if (body) {
		SET_VEC(v,rp->spin[rp->mom_ord[X]],rp->spin[rp->mom_ord[Y]],rp->spin[rp->mom_ord[Z]]);
		COPY_MAT(rp->axes,m);
		Transpose(m);
		Transform(m,v,rp->spin);
		}

	for (i=0;i<rp->n_particles;i++) {
		d = &rp->data[i];
		SUB_VEC(d->pos,rp->pos,r);
		CROSS(rp->spin,r,v);
		ADD_VEC(d->vel,v,d->vel);
		ADD_VEC(d->spin,rp->spin,d->spin);
		}
	}

void rpuAddAngMom(RUBBLE_PILE *rp)
{
	/*
	 ** Augments particle velocities by specified angular momentum
	 ** by computing equivalent spin and calling rpuAddSpin(). Note the
	 ** inertia tensor is required in advance -- see rpuCalcInertia().
	 */

	MATRIX a;

	assert(rp != NULL);
	assert(rp->n_particles > 0);
	assert(rp->data != NULL);

	assert(rp->mass >= 0.0);

	if (rp->mass > 0.0) {
		COPY_MAT(rp->inertia,a);
		invert(a);
		Transform(a,rp->ang_mom,rp->spin);
		}

	rpuAddSpin(rp,FALSE);
	}

void rpuAnalyze(RUBBLE_PILE *rp)
{
	/*
	 ** Given rp->n_particles and rp->data, this routine fills in the
	 ** remaining rubble pile parameters, which are basically global
	 ** properties. Note the calculation order is important.
	 */

	rpuCalcMass(rp);
	rpuCalcPos(rp);
	rpuCalcVel(rp);
	rpuCalcRadius(rp);
	rpuCalcInertia(rp);
	rpuCalcAxes(rp);
	rpuCalcSpin(rp);
	rpuCalcColor(rp);
	rpuCalcAggID(rp);
	rpuCalcDensity(rp);
	rpuCalcPacking(rp);
	}

/*
 ** Following routines adapted from Numerical Recipes (3rd ed).
 */

static void eigsrt(VECTOR d,MATRIX v)
{
	/* based on eigsrt(), NR3 11.1 */

	const int n = N_DIM; /* i.e. hard-wired for 3-D */

	int i,j,k;
	double p;

	for (i=0;i<n-1;i++) {
		p=d[k=i];
		for (j=i;j<n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=0;j<n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
				}
			}
		}
	}

/*inline*/ void rot(MATRIX a,double s,double tau,int i,int j,int k,int l)
{
	/* from struct Jacobi, NR3 11.1 */

	double g=a[i][j];
	double h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
	}

static void jacobi(MATRIX a,VECTOR d,MATRIX v)
{
	/* based on struct jacobi(), NR3 11.1 */

	const int n = N_DIM; /* i.e. hard-wired for 3-D */
	const double EPS = 1.1102e-16; /* numeric_limits<Doub>::epsilon() */

	int i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;
	VECTOR b,z;

	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
		}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
		}
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
			}
		if (sm == 0.0) {
			eigsrt(d,v);
			return;
			}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && g <= EPS*fabs(d[ip]) && g <= EPS*fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (g <= EPS*fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
						}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;j++)
						rot(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;j++)
						rot(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;j++)
						rot(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;j++)
						rot(v,s,tau,j,ip,j,iq);
					}
				}
			}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
			}
		}
	(void) fprintf(stderr, "jacobi(): Too many iterations.\n");
	exit(1);
	}

/*
 ** Following routines adapted from Numerical Recipes in C (2nd ed).
 */

static void ludcmp(MATRIX a,int *indx)
{
	/* based on ludcmp(), NRiC(2e) 2.3 */

	const int n = N_DIM;

	VECTOR vv;
	double big,dum,sum,temp;
	int imax=0,i,j,k;

	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) {
			(void) fprintf(stderr, "ludcmp(): Singular matrix.\n");
			exit(1);
			}
		vv[i]=1.0/big;
		}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
				}
			}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
				}
			vv[imax]=vv[j];
			}
		indx[j]=imax;
		/*if (a[j][j] == 0.0) a[j][j]=TINY;*/
		if (j != n - 1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
			}
		}
	}

static void lubksb(MATRIX a,int *indx,VECTOR b)
{
	/* based on lubksb(), NRiC(2e) 2.3 */

	const int n = N_DIM;

	double sum;
	int ip,ii=(-1),i,j;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii >= 0)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
		}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
		}
	}

static void invert(MATRIX a)
{
	/* inverts matrix */

	MATRIX b;
	VECTOR col;
	int indx[N_DIM],i,j;

	ludcmp(a,indx);

	for (j=0;j<N_DIM;j++) {
		for (i=0;i<N_DIM;i++)
			col[i] = 0;
		col[j] = 1;
		lubksb(a,indx,col);
		for (i=0;i<N_DIM;i++)
			b[i][j] = col[i];
		}

	COPY_MAT(b,a);
	}

/* rpu.c */

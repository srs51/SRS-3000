/*
 ** ssgen.c -- DCR 02/03/12
 ** =====
 ** Generates initial conditions in a rectangular or ellipsoidal
 ** space, allowing for uncertainty in particle size, and using a tree
 ** code to eliminate particle overlaps.
 ** update -- RLB 06/04/12
 ** now includes bi-modal and tri-modal radius and density distribution, and
 ** discontinuous and continuous power-law distribution 
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h> /* for getpid() */
#include <unistd.h> /* for getopt() and getpid() */
#include <string.h> /* for strlen() */
#include <ctype.h> /* for tolower() */
#include <time.h>
#include <math.h>
#include <assert.h>
#include <ss.h>
#include <random.h>

#define OUTFILENAME "ssgen.ss"
#define LOGFILENAME "ssgen.log"

#define MAX_NUM_PASS 500

typedef int BOOLEAN;

typedef struct {                
  int nPart,bEllipsoidal;  
  double vAxes[N_DIM],dRadAvg,dBRad,dBDensity,dTRad,dTDensity,dRadDev,dCutoff,dDensity,dPIndex,dPRmax,dPFlag;
  } PARAMS;

// position of particle, and particle radius? would need to ammend this for bimodality etc... -rb
typedef struct {
  int keep;
  double pos[N_DIM],radius,density;
  } DATA;

double SQ(double x);
#define SQ(x) ((x)*(x))

double CUBE(double x);
#define CUBE(x) ((x)*SQ(x))

double DISTSQ(double v1[3],double v2[3]);
#define DISTSQ(v1,v2) ((SQ(v1[0] - v2[0]) + SQ(v1[1] - v2[1]) + SQ(v1[2] - v2[2])))

BOOLEAN get_yn(const char *str,const char *dflt_str)
{
	enum {none,yes,no} dflt = none;

	char c;

	if (dflt_str && strlen(dflt_str)) {
		if (tolower(*dflt_str) == 'y') dflt = yes;
		else if (tolower(*dflt_str) == 'n') dflt = no;
		}

	do {
		(void) printf("%s [%s]? ",str,
					  (dflt == none ? "y/n" : dflt == yes ? "Y/n" : "y/N"));
		c = tolower(getchar());
		if (c == '\n' && dflt != none)
			return dflt == yes;
		else if (c == 'y' || c == 'n') {
			BOOLEAN is_yes = (c == 'y');
			do c = getchar(); /* eat any leftover characters */
			while (c != '\n');
			return is_yes;
			}
		while (c != '\n') c = getchar();
		} while (/*CONSTCOND*/1);
	}

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  double cmin[N_DIM],cmin_eff[N_DIM],cmax[N_DIM],cmax_eff[N_DIM],pos[N_DIM]; //NDIM=3 defined in ssio.h
  struct node_s *cell[CELLS_PER_NODE];
  DATA *leaf[CELLS_PER_NODE];
  };

typedef struct node_s NODE;

//writing the logfile - rb
static void write_log(const PARAMS *p,const DATA *d)
{
	FILE *fp;

	fp = fopen(LOGFILENAME,"w");
	if (fp == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",LOGFILENAME);
		exit(1);
		}

	fprintf(fp,"Target particle radius = %g cm\n",100.0*p->dRadAvg*AU);
	// add if statement here if bi-modality or tri-modality is specified.
	if (p->dBRad > 1 && p->dTRad == 1.0){	
		fprintf(fp,"Bi-modal Distribution specified with:\n");
		fprintf(fp,"	Radius Ratio = %.2f, Density Ratio = %.2f \n",p->dBRad,p->dBDensity);}
	if (p->dBRad > 1 && p->dTRad > 1){	
		fprintf(fp,"Tri-modal Distribution specified with:\n");
		fprintf(fp,"	Radius Ratio = %.2f:%.2f:1, Density Ratio = %.2f:%.2f:1 \n",p->dTRad,p->dBRad,p->dTDensity,p->dBDensity);}
	if (p->dPFlag == 1){
		fprintf(fp,"Discontinuous Power-Law Distribution specified with:\n");
		fprintf(fp,"	Index = %.2f, Rmax/Rmin = %.2f \n",p->dPIndex,p->dPRmax);}
	if (p->dPFlag == 2){
		fprintf(fp,"Continuous Power-Law Distribution specified with:\n");
		fprintf(fp,"	Index = %.2f, Rmax/Rmin = %.2f \n",p->dPIndex,p->dPRmax);}
	fprintf(fp,"Index of power law distribution = %g \n",p->dPIndex);
	fprintf(fp,"Target particle radius uncertainty = %g cm\n",100.0*p->dRadDev*AU);
	fprintf(fp,"Gaussian distribution cutoff = %g sigma\n",p->dCutoff);
	fprintf(fp,"Particle density = %g g/cc\n",1.0e-3*p->dDensity*M_SUN/(AU*AU*AU));
	fprintf(fp,"Region dimensions = %g x %g x %g m (%s)\n",p->vAxes[0]*AU,p->vAxes[1]*AU,p->vAxes[2]*AU,p->bEllipsoidal ? "ellipsoidal" : "rectangular");
	fprintf(fp,"Number of particles = %i\n",p->nPart);
	//probably need to change this as this now reflects the average radius of the smallest particle
	if (p->nPart > 1) {
		double x,s;
		int i;

		x = 0.0;
		for (i=0;i<p->nPart;i++) // determining the average particle radius by adding up all the radii and diving by the # of particles
			x += d[i].radius;
		x /= p->nPart;
		s = 0.0;
		for (i=0;i<p->nPart;i++) // determining the variance in particle radius 
			s += SQ(d[i].radius - x); // add up all the differences between each particle radius and average 
		s = sqrt(s/(p->nPart - 1)); // take sqrt of (value / (# of particles -1))
		fprintf(fp,"Average particle radius = %g +/- %g cm\n",100.0*x*AU,100.0*s*AU);
		}

	fprintf(fp,"Minimum particle radius = %g cm\n",100.0*d[0].radius*AU); //particles indices are arranged by radius, determined by other function?
	fprintf(fp,"Maximum particle radius = %g cm\n",100.0*d[p->nPart - 1].radius*AU);
	fprintf(fp,"Suggested minimum nSmooth:\n");
	/* for HSDEM, ratio of largest particle surface area to smallest particle cross-section */
	fprintf(fp,"   HSDEM = %i\n",4*((int) (SQ(d[p->nPart - 1].radius/d[0].radius) + 0.5)));	//nsmooth  = 4*(sqrt(largest_radius/smallest_radius + 1/2)) - an integer
	/* for SSDEM, number of small particles in sphere of radius small + large particle radii */
	fprintf(fp,"   SSDEM = %i\n",(int) (CUBE(d[p->nPart - 1].radius/d[0].radius + 1.0) + 0.5));
	fprintf(fp,"   (note: a smaller value for SSDEM is possible depending on filling factor)\n");
	fclose(fp);
	}
//function that is responsible for writing to ssgen.ss 
static void write_data(const PARAMS *p,const DATA *d)
{
	// these data-structures are defined in ssio.h 
	SSIO ssio;	
	SSHEAD head;
	SSDATA data;
	int i;
	//checks that ssgen.ss is open for writing?
	if (ssioOpen(OUTFILENAME,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",OUTFILENAME);
		exit(1);
		}

	head.time = 0.0;
	head.n_data = p->nPart;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;
	//makes sure header is properly written (ss standard header?) ?
	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Unable to write header.\n");
		ssioClose(&ssio);
		exit(1);
		}
	//writing data... to ssdata data structure
	//using top defined DATA structure to define SSDATA members (named data.), defined within this function
	for (i=0;i<p->nPart;i++) {
		data.mass = 4.0/3.0*M_PI*d[i].radius*d[i].radius*d[i].radius*d[i].density; //mass given by particle radius and density (given by user)
		data.radius = d[i].radius;
		data.pos[0] = d[i].pos[0];
		data.pos[1] = d[i].pos[1];
		data.pos[2] = d[i].pos[2];
		data.vel[0] = data.vel[1] = data.vel[2] = 0.0;
		data.spin[0] = data.spin[1] = data.spin[2] = 0.0;
		data.color = PLANETESIMAL;
		data.org_idx = i;
		if (ssioData(&ssio,&data)) {
			fprintf(stderr,"Error writing data (particle %i).\n",i);
			ssioClose(&ssio);
			exit(1);
			}
		}

	ssioClose(&ssio);
	}
//node's used for rejection algorithm 
static void kill_node(NODE *node)
{
	int i;

	assert(node != NULL);

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i])
			kill_node(node->cell[i]);

	free((void *) node);
	}
//function that squezzes particles into box by checking for overlaps. Uses node_s structure
static int reject(const NODE *node,const DATA *d)
{
	NODE *n;
	double r2;
	int i,k;

	/* intersect test -- check each subcell in turn */

	/*
	** Algorithm adapted from...
	** "A Simple Method for Box-Sphere Intersection Testing",
	** by Jim Arvo, in "Graphics Gems", Academic Press, 1990.
	** (http://www.ics.uci.edu/~arvo/code/BoxSphereIntersect.c)
	*/

	for (i=0;i<CELLS_PER_NODE;i++) {
		n = node->cell[i];
		if (n != NULL) {
			r2 = 0.0;
			for (k=0;k<N_DIM;k++)
				if (d->pos[k] < n->cmin_eff[k])
					r2 += SQ(d->pos[k] - n->cmin_eff[k]);
				else if (d->pos[k] > n->cmax_eff[k])
					r2 += SQ(d->pos[k] - n->cmax_eff[k]);
			if (r2 <= SQ(d->radius) && reject(n,d)) /*reject called through root than loops through all the nodes till it finds a NULL*/
				return 1;
			}
		else if (node->leaf[i] != NULL && node->leaf[i] != d &&
				 node->leaf[i]->keep &&
				 DISTSQ(d->pos,node->leaf[i]->pos) <=
				 SQ(d->radius + node->leaf[i]->radius))
			return 1;
		}

	return 0;
	}
// makes node given xyz min's and maxes from user input, to be used by rejection algorithm
static void make_node(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,NODE **node)
{
	int i;

	*node = (NODE *) malloc(sizeof(NODE));
	assert(*node != NULL);

	(*node)->cmin[0] = (*node)->cmin_eff[0] = xmin; /* define boundaries of rectangular cell */
	(*node)->cmin[1] = (*node)->cmin_eff[1] = ymin;
	(*node)->cmin[2] = (*node)->cmin_eff[2] = zmin;
	(*node)->cmax[0] = (*node)->cmax_eff[0] = xmax;
	(*node)->cmax[1] = (*node)->cmax_eff[1] = ymax;
	(*node)->cmax[2] = (*node)->cmax_eff[2] = zmax;

	(*node)->pos[0] = 0.5*((*node)->cmin[0] + (*node)->cmax[0]); /* define center of cell */
	(*node)->pos[1] = 0.5*((*node)->cmin[1] + (*node)->cmax[1]);
	(*node)->pos[2] = 0.5*((*node)->cmin[2] + (*node)->cmax[2]);

	for (i=0;i<CELLS_PER_NODE;i++) { /* initalize the components of the struct, also required to end recursion */
		(*node)->cell[i] = NULL;
		(*node)->leaf[i] = NULL;
		}
	}
// part of rejection algorithm
static void add_to_tree(NODE *node,DATA *d)
{
	int i,idx,idy,idz;

	idx = (d->pos[0] < node->pos[0] ? -1 : 1); //tri-conditional
	idy = (d->pos[1] < node->pos[1] ? -1 : 1);
	idz = (d->pos[2] < node->pos[2] ? -1 : 1);

	i = (idx + 1)/2 + (idy + 1 + 2*(idz + 1)); /* determine which octant the node/particle belongs in */

	if (node->cell[i] != NULL)
		add_to_tree(node->cell[i],d);
	else if (node->leaf[i] != NULL) {
		make_node(idx < 0 ? node->cmin[0] : node->pos[0],
				  idx < 0 ? node->pos[0] : node->cmax[0],
				  idy < 0 ? node->cmin[1] : node->pos[1],
				  idy < 0 ? node->pos[1] : node->cmax[1],
				  idz < 0 ? node->cmin[2] : node->pos[2],
				  idz < 0 ? node->pos[2] : node->cmax[2],
				  &node->cell[i]);
		add_to_tree(node->cell[i],node->leaf[i]);
		add_to_tree(node->cell[i],d);
		node->leaf[i] = NULL;
		}
	else {
		node->leaf[i] = d;
		if (node->cmin[0] - d->radius < node->cmin_eff[0])
			node->cmin_eff[0] = node->cmin[0] - d->radius;
		if (node->cmax[0] + d->radius > node->cmax_eff[0])			
			node->cmax_eff[0] = node->cmax[0] + d->radius;
		if (node->cmin[1] - d->radius < node->cmin_eff[1])
			node->cmin_eff[1] = node->cmin[1] - d->radius;
		if (node->cmax[1] + d->radius > node->cmax_eff[1])
			node->cmax_eff[1] = node->cmax[1] + d->radius;
		if (node->cmin[2] - d->radius < node->cmin_eff[2])
			node->cmin_eff[2] = node->cmin[2] - d->radius;
		if (node->cmax[2] + d->radius > node->cmax_eff[2])
			node->cmax_eff[2] = node->cmax[2] + d->radius;
			
		}
	}
//making sure that radius matchup in data. (ordering by radius?)
int compar(const void *v1,const void *v2)
{
	if (((DATA *)v1)->radius < ((DATA *)v2)->radius)
		return -1;
	else if (((DATA *)v1)->radius > ((DATA *)v2)->radius)
		return 1;
	else
		return 0;
	}
/*gives random position to particle. Does so once if rectangular, loops through if ellipsoidal and not within ellipse constraints*/
void get_pos(const PARAMS *p,DATA *d)
{
	/*
	** NOTE: for the ellipsoidal region, this routine guarantees that
	** the particle center -- but not necessarily the entire particle
	** -- is contained within the ellipsoid (both are guaranteed for
	** the rectangular region).
	*/

	int k;
	double a;

	do {
		for (k=0;k<N_DIM;k++) {
			a = p->vAxes[k] - 2.0*d->radius;
			if (a < 0.0) {
				fprintf(stderr,"get_pos(): Region too small for particle.\n");
				exit(1);
				}
			d->pos[k] = (randUniform() - 0.5)*a;
			}			
		} while (p->bEllipsoidal &&
				 SQ(d->pos[0])/SQ(p->vAxes[0]) +
				 SQ(d->pos[1])/SQ(p->vAxes[1]) +
				 SQ(d->pos[2])/SQ(p->vAxes[2]) > 0.25);
	}

double gaussian(double dSig,double dLim)
{
	/* truncated Gaussian */

	double dVal;

	do {
		dVal = dSig*randGaussian();
		} while (dLim > 0.0 && fabs(dVal) > dLim*dSig);

	return dVal;
	}
/* generating particles */
static void generate(const PARAMS *p,DATA *d)
{
	NODE *root;
	BOOLEAN bTest;
	long i,j,k=0,res;
	int nPass,nRej;
	double Nb,Nm,Ns,unidev;
	randSeedGenerator(time(NULL) % getpid() + getppid());

	/* pregenerate particle radii */
	res=(long)p->nPart; /* set resolution to total number of particles (i.e. each particle has a unique radius)
			      Seems to work best for smooth distribution. If this is optimal, then remove second for loop
			      Currently keeping this here if a resolution change is desired in future (possibly as user-input)
			      However, may be best to remove resolution parameter completely so as to remove possibility of resolution > N_particles*/
	/*Continuous and Discontinuous Power Law Distribution*/
	if (p->dPFlag == 2.0){
				
		if (p->dPIndex != -1){	
		for (i=0;i<res;i++){
			for (j=k;j<p->nPart*(i+1)/res;j++){
				d[j].keep=0;
			do{
				unidev=(float) (i+1)/res;
				d[j].radius=pow(((1-unidev)*pow(p->dRadAvg,p->dPIndex +1))+(unidev*pow(p->dPRmax*p->dRadAvg,p->dPIndex+1)),1/(p->dPIndex+1));
				d[j].density=p->dDensity;
				} while (d[j].radius <= 0.0);			
			}
			k=j;
		}}

		if (p->dPIndex == -1){
		for (i=0;i<res;i++){
			for (j=k;j<p->nPart*(i+1)/res;j++){
				d[j].keep=0;
			do{
				unidev=(float) (i+1)/res;
				d[j].radius=pow(p->dPRmax*p->dRadAvg,unidev)/pow(p->dRadAvg,unidev-1);
				d[j].density=p->dDensity;
				} while (d[j].radius <= 0.0);			
			}
			k=j;
		}}				
	}
	
	if (p->dPFlag == 1.0){
		if (p->dPIndex != -1){	
		for (i=0;i<p->nPart;i++){
			d[i].keep = 0;
			do{	
				unidev=randUniform();
				d[i].radius=pow(((1-unidev)*pow(p->dRadAvg,p->dPIndex +1))+(unidev*pow(p->dPRmax*p->dRadAvg,p->dPIndex+1)),1/(p->dPIndex+1));
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}}

		if (p->dPIndex == -1){
		for (i=0;i<p->nPart;i++){
			d[i].keep = 0;
			do{	
				unidev=randUniform();
				d[i].radius=pow(p->dPRmax*p->dRadAvg,unidev)/pow(p->dRadAvg,unidev-1);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}}

	}

	//Whether or not bi-modality is slected this will generate particle radii based on equal Volume
	if (p->dBRad > 1.0){
		bTest = get_yn("Bimodality: Specify Number of Larger Particles?","n");
		if (bTest){	
			do {
				printf("Enter Number of Larger Particles: ");
				(void) scanf("%lf",&Nm);
				if (Nm == 0.0) printf("Number of Larger Particles cannot be zero.\n");
				} while (Nm == 0.0);
		}	
		else Nm=p->nPart/(1+pow(p->dBRad,3));
		

		Ns=p->nPart-Nm;

		for (i=0;i<Ns;i++) {
			d[i].keep = 0;
			do {
				d[i].radius = p->dRadAvg + gaussian(p->dRadDev,p->dCutoff);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}
		for (j=i;j<p->nPart;j++) {
			d[j].keep = 0;
			do {
				d[j].radius = (p->dBRad)*(p->dRadAvg + gaussian(p->dRadDev,p->dCutoff));
				d[j].density=(p->dBDensity)*p->dDensity;
				} while (d[j].radius <= 0.0); /* reject unphysical values */
			}
	}

	if (p->dTRad > 1.0){
		Nb=p->nPart/( 1+pow(p->dTRad,3)+pow(p->dTRad/p->dBRad,3));
		Nm=p->nPart/( 1+pow(p->dBRad,3)+pow(p->dBRad/p->dTRad,3));
		Ns=p->nPart-Nm-Nb;

		for (i=0;i<Ns;i++) {
			d[i].keep = 0;
			do {
				d[i].radius = p->dRadAvg + gaussian(p->dRadDev,p->dCutoff);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}
		for (j=i;j<Ns+Nm;j++) {
			d[j].keep = 0;
			do {
				d[j].radius = (p->dBRad)*(p->dRadAvg + gaussian(p->dRadDev,p->dCutoff));
				d[j].density=(p->dBDensity)*p->dDensity;
				} while (d[j].radius <= 0.0); /* reject unphysical values */
			}

		for (k=j;k<p->nPart;k++) {
			d[k].keep = 0;
			do {
				d[k].radius = (p->dTRad)*(p->dRadAvg + gaussian(p->dRadDev,p->dCutoff));
				d[k].density=(p->dTDensity)*p->dDensity;
				} while (d[k].radius <= 0.0); /* reject unphysical values */
			}
	}

	//default arrangment
	if (p->dTRad == 1.0 && p->dBRad == 1.0 && p->dPFlag == 0.0){
		for (i=0;i<p->nPart;i++) {
			d[i].keep = 0;
			do {
				d[i].radius = p->dRadAvg + gaussian(p->dRadDev,p->dCutoff);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}
	}
		
	//putting this here as a precaution for now
        for (i=0;i<p->nPart;i++) d[i].keep=0;
		

	/* sort particles by increasing radius */
	qsort(d,p->nPart,sizeof(DATA),compar);
	
	// rejection algorithm
	nPass = 0;
	do {
		if (++nPass > MAX_NUM_PASS) {
			fprintf(stderr,"Unable to converge after %i pass%s.\n",
					nPass - 1,nPass == 2 ? "" : "es");
			exit(1);
			}
		printf("Pass %i ",nPass);
		if (nPass == 1)
			printf("(nPart = %i)\n",p->nPart);
		else
			printf("(nRej = %i)\n",nRej);
		/* make root cell */
		make_node(-0.5*p->vAxes[0],0.5*p->vAxes[0],
				  -0.5*p->vAxes[1],0.5*p->vAxes[1],
				  -0.5*p->vAxes[2],0.5*p->vAxes[2],&root);
		/* generate particle positions (as needed) and add to tree */
		for (i=0;i<p->nPart;i++) {
			if (!d[i].keep) {
				get_pos(p,&d[i]);
				d[i].keep = 1;
				}
			add_to_tree(root,&d[i]);
			}
		/* check for overlaps, rejecting smallest particles first */
		nRej = 0;
		for (i=0;i<p->nPart;i++)
			if (reject(root,&d[i])) {
				d[i].keep = 0;
				++nRej;
				}
		kill_node(root);
		} while (nRej > 0);
	}

static void usage(const char *achProgName)
{
	fprintf(stderr,
			"Usage: %s -r radius -d density [-b bi-modal RadRatio,DensityRatio] [ -t tri-modal R3:R2,D3:D2] [-p power-law index,Rmax/Rmin,flag] [ -s stddev ] [ -q cutoff ]  -x lx -y ly -z lz [ -e ] N\n"
			"where: radius is the mean particle radius (cm) [min radius for -b, -t, or -p];\n"
			"       bi-modal radius (>= 1.0) &/or density distribution (> 0), e.g. -b 1.4,0.8;\n"
			"       tri-modal distribution, e.g. -t 2.4:1.2,1.1:3.5 gives:\n"
			"          3 types of particles...one has 2.4x radius & 1.1x density,\n"
			"          2nd has 1.2x radius & 3.5x density,\n"
			"          3rd has radius of -r and density of -d;\n"
			"       power-law radius: power law index (-1 included), Max to Min radius,\n"
			"          flag for discontinuous (=1) or continuous (=2) distribution;\n"
			"       stddev is the uncertainty in particle radius (cm);\n"
			"       cutoff is max. no. sigma for dispersions (default 1; 0 means no limit);\n"
			"       density is the particle density (g/cc);\n"
			"       lx, ly, lz are the simulation region dimensions (m);\n"
			"       -e means ellisoidal, not rectangular, region;\n"
			"       and N is the number of particles to generate.\n"
			"       NOTE: region dimensions are full axes, not semi axes;\n"
			"             particle radius taken into account where possible.\n",
			achProgName);

	exit(1);
}

int main(int argc,char *argv[])
{
	extern char *optarg; //user input
	extern int optind; //current index into argument array
	
	PARAMS params;
	DATA *data;
	double d,dd,ddd,dddd;
	char *rstring, *dstring, *r3, *r2, *d3, *d2, *istring, *fstring;	
	int c;
	
	assert(MAX_NUM_PASS > 0);

	/* defaults */

	params.dBRad = 1.0; //new parameter need to define
	params.dBDensity = 0.0;
	params.dTRad = 1.0;
	params.dTDensity = 0.0;
	params.dPIndex = 0.0;
	params.dPRmax = 1.0;
	params.dPFlag = 0.0;

	params.dRadAvg = 0.0;
	params.dRadDev = 0.0;
	params.dCutoff = 1.0;
	params.dDensity = 0.0;
	params.vAxes[0] = params.vAxes[1] = params.vAxes[2] = 0.0;
	params.bEllipsoidal = 0;
	params.nPart = 0;

	/* parse command-line arguments */

	while ((c = getopt(argc,argv,"r:b:t:p:s:q:d:x:y:z:e")) != EOF)
		switch (c) {
		case 'r':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dRadAvg = 0.01*d/AU; /* cm -> AU */
			break;
		case 'b':
			rstring = strtok (optarg, ","); //parsing input string
			dstring = strtok (NULL, ",");
			d = atof(rstring);
			dd = atof(dstring);
			if (d < 1.0)
				usage(argv[0]);
			if (dd <= 0.0)
				usage(argv[0]);
			
			params.dBRad = d; 
			params.dBDensity =dd; 
			break;
		case 't':
			r3 = strtok (optarg, ":");
			r2 = strtok (NULL, ",");
			d3 = strtok (NULL, ":");
			d2= strtok (NULL, ":");
			d = atof(r2);
			dd = atof(d2);
			ddd = atof(r3);
			dddd = atof (d3);
			if (d <= 1.0)
				usage(argv[0]);
			if (dd <= 0.0)
				usage(argv[0]);
			if (ddd <= 1.0)
				usage(argv[0]);
			if (dddd <= 0.0)
				usage(argv[0]);
			
			params.dBRad = d; 
			params.dBDensity = dd;
			params.dTRad = ddd;
			params.dTDensity = dddd; 
			break;
		case 'p':
			istring = strtok (optarg, ",");
			rstring = strtok (NULL, ",");
			fstring = strtok (NULL, ",");
			d = atof(istring);
			dd = atof(rstring);
			ddd = atof(fstring); 			
			if (dd <= 1.0)
				usage(argv[0]);
			if (ddd != 1 && ddd != 2)
				usage(argv[0]);
			params.dPIndex = d; 
			params.dPRmax = dd;
			params.dPFlag = ddd; 
			break;				
		case 's':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dRadDev = 0.01*d/AU; /* cm -> AU */
			break;
		case 'q':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dCutoff = d;
			break;
		case 'd':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dDensity = 1.0e3*d/M_SUN*AU*AU*AU; /* g/cc -> M_SUN/AU^3 */
			break;
		case 'x':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.vAxes[0] = d/AU; /* m -> AU */
			break;
		case 'y':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.vAxes[1] = d/AU; /* m -> AU */
			break;
		case 'z':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.vAxes[2] = d/AU; /* m -> AU */
			break;
		case 'e':
			params.bEllipsoidal = 1;
			break;
		case '?':
		default:
			usage(argv[0]);
			}

	if (optind >= argc || params.dRadAvg == 0.0 || params.dDensity == 0.0 || params.vAxes[0] == 0.0 || params.vAxes[1] == 0.0 || params.vAxes[2] == 0.0)
		usage(argv[0]);

	if (params.dRadDev >= params.dRadAvg)
		fprintf(stderr,"WARNING: large particle radius uncertainty---infinite loop possible.\n");

	params.nPart = atoi(argv[optind]);

	if (params.nPart <= 0 || argc > optind + 1)
		usage(argv[0]);

	data = (DATA *) malloc(params.nPart*sizeof(DATA));
	assert(data != NULL);

	generate(&params,data);

	write_log(&params,data);
	printf("Log written to %s.\n",LOGFILENAME);

	write_data(&params,data);

	free(data);

	printf("Now run rpx on %s to reset center of mass, etc.\n",OUTFILENAME);

	return 0;
	}

/* ssgen.c */

/*
 ** ssn.c -- DCR 01/23/12
 ** =====
 ** Outputs neighbor lists of ss particles.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() */
#include <math.h>
#include <assert.h>
#include <ssio.h>

#define OUT_EXT ".nbr"

typedef struct {
  double dSearchRadius;
  double dSearchRadius2;
  } PARAMS;

typedef struct {
  int iIndex;
  double dMass,pos[N_DIM];
  } DATA;

double SQ(double x);
#define SQ(x) ((x)*(x))

double DISTSQ(double v1[3],double v2[3]);
#define DISTSQ(v1,v2) ((SQ(v1[0] - v2[0]) + SQ(v1[1] - v2[1]) + SQ(v1[2] - v2[2])))

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  double cmin[N_DIM],cmax[N_DIM],pos[N_DIM];
  struct node_s *cell[CELLS_PER_NODE];
  DATA *leaf[CELLS_PER_NODE];
  };

typedef struct node_s NODE;

static void kill_node(NODE *node)
{
	int i;

	assert(node != NULL);

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i])
			kill_node(node->cell[i]);

	free((void *) node);
	}

static void search(const PARAMS *p,const DATA *d,const NODE *node,FILE *fp)
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
				if (d->pos[k] < n->cmin[k])
					r2 += SQ(d->pos[k] - n->cmin[k]);
				else if (d->pos[k] > n->cmax[k])
					r2 += SQ(d->pos[k] - n->cmax[k]);
			if (r2 <= p->dSearchRadius2)
				search(p,d,n,fp);
			}
		else if (node->leaf[i] != NULL && node->leaf[i] != d &&
				 (r2 = DISTSQ(d->pos,node->leaf[i]->pos)) <= p->dSearchRadius2)
			fprintf(fp," %i %e",node->leaf[i]->iIndex,p->dSearchRadius - sqrt(r2));
		}
	}

static void make_node(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,NODE **node)
{
	int i;

	*node = (NODE *) malloc(sizeof(NODE));
	assert(*node != NULL);

	(*node)->cmin[0] = xmin;
	(*node)->cmax[0] = xmax;
	(*node)->cmin[1] = ymin;
	(*node)->cmax[1] = ymax;
	(*node)->cmin[2] = zmin;
	(*node)->cmax[2] = zmax;

	(*node)->pos[0] = 0.5*((*node)->cmin[0] + (*node)->cmax[0]);
	(*node)->pos[1] = 0.5*((*node)->cmin[1] + (*node)->cmax[1]);
	(*node)->pos[2] = 0.5*((*node)->cmin[2] + (*node)->cmax[2]);

	for (i=0;i<CELLS_PER_NODE;i++) {
		(*node)->cell[i] = NULL;
		(*node)->leaf[i] = NULL;
		}
	}

static void add_to_tree(NODE *node,DATA *d)
{
	int i,idx,idy,idz;

	idx = (d->pos[0] < node->pos[0] ? -1 : 1);
	idy = (d->pos[1] < node->pos[1] ? -1 : 1);
	idz = (d->pos[2] < node->pos[2] ? -1 : 1);

	i = (idx + 1)/2 + (idy + 1 + 2*(idz + 1));

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
	else
		node->leaf[i] = d;
	}

static void find_nbrs(const PARAMS *p,DATA d[],int nData,FILE *fp)
{
	NODE *root;
	double xcom,ycom,zcom,xmax,ymax,zmax;
	double m,dMass;
	int i;

	/* get center of mass */

	xcom = ycom = zcom = dMass = 0.0;
	for (i=0;i<nData;i++) {
		m = d[i].dMass;
		xcom += m*d[i].pos[0];
		ycom += m*d[i].pos[1];
		zcom += m*d[i].pos[2];
		dMass += m;
	}

	assert(dMass > 0.0);
	xcom /= dMass;
	ycom /= dMass;
	zcom /= dMass;

	/* adjust center of mass to origin and get root cell dimensions */

	xmax = ymax = zmax = -1.0;
	for (i=0;i<nData;i++) {
		d[i].pos[0] -= xcom;
		d[i].pos[1] -= ycom;
		d[i].pos[2] -= zcom;
		m = fabs(d[i].pos[0]);
		if (m > xmax)
			xmax = m;
		m = fabs(d[i].pos[1]);
		if (m > ymax)
			ymax = m;
		m = fabs(d[i].pos[2]);
		if (m > zmax)
			zmax = m;
		}

	/* allocate and initialize root cell */

	make_node(-xmax,xmax,-ymax,ymax,-zmax,zmax,&root);

	/* fill tree with particles */

	for (i=0;i<nData;i++)
		add_to_tree(root,&d[i]);

	/* find neighbors of each particle */

	fprintf(fp,"%i\n",nData);
	for (i=0;i<nData;i++) {
		fprintf(fp,"%i",i);
		search(p,&d[i],root,fp);
		fprintf(fp,"\n");
		}

	/* deallocate tree */

	kill_node(root);
	}

static int read_data(const char *achFile,DATA **d,int *nData)
{
	SSIO ssio;
	SSHEAD h;
	SSDATA ssd;
	int i;

	printf("%s: ",achFile);

	if (ssioOpen(achFile,&ssio,SSIO_READ)) {
		fprintf(stderr,"Unable to open \"%s\" for reading.\n",achFile);
		return 1;
		}

	if (ssioHead(&ssio,&h)) {
		fprintf(stderr,"Corrupt header.\n");
		ssioClose(&ssio);
		return 1;
		}

	if (h.n_data <= 0) {
		fprintf(stderr,"Invalid file size.\n");
		ssioClose(&ssio);
		return 1;
		}

	if (h.iMagicNumber != SSIO_MAGIC_STANDARD &&
		h.iMagicNumber != SSIO_MAGIC_REDUCED) {
		fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio);
		return 1;
		}
	
	*nData = h.n_data;
	*d = (DATA *) malloc(*nData*sizeof(DATA));
	if (*d == NULL) {
		fprintf(stderr,"Unable to allocate memory for %i particle%s.\n",*nData,
				*nData == 1 ? "" : "s");
		return 1;
		}

	for (i=0;i<*nData;i++) {
		if (ssioData(&ssio,&ssd)) {
			fprintf(stderr,"Corrupt data\n");
			free(*d);
			(void) ssioClose(&ssio);
			return 1;
			}
		(*d)[i].iIndex = i;
		(*d)[i].dMass = ssd.mass;
		(*d)[i].pos[0] = ssd.pos[0];
		(*d)[i].pos[1] = ssd.pos[1];
		(*d)[i].pos[2] = ssd.pos[2];
		}

	printf("read %i particle%s.\n",*nData,*nData == 1 ? "" : "s");

	ssioClose(&ssio);

	return 0;
	}

static void usage(const char *achProgName)
{
	fprintf(stderr,
			"Usage: %s -R radius ss-file [ ss-file ... ]\n"
			"where radius is the search radius to use (pkdgrav units).\n",
			achProgName);
	exit(1);
}

int main(int argc,char *argv[])
{
	extern char *optarg;
	extern int optind;

	PARAMS params;
	DATA *data;
	FILE *fp;
	char achOutFile[MAXPATHLEN];
	double d;
	int c,nData;

	params.dSearchRadius2 = 0.0;

	/* parse command-line arguments */

	while ((c = getopt(argc,argv,"R:")) != EOF)
		switch (c) {
		case 'R':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dSearchRadius = d;
			params.dSearchRadius2 = SQ(d);
			break;
		case '?':
		default:
			usage(argv[0]);
			}

	if (optind >= argc || params.dSearchRadius2 == 0.0)
		usage(argv[0]);

	/* loop over ss files */

	for (;optind<argc;optind++) {
		if (read_data(argv[optind],&data,&nData) != 0)
			continue;
		ssioNewExt(argv[optind],SS_EXT,achOutFile,OUT_EXT);
		fp = fopen(achOutFile,"w");
		if (fp == NULL) {
			fprintf(stderr,"Unable to open \"%s\" for writing.\n",achOutFile);
			free(data);
			continue;
			}
		find_nbrs(&params,data,nData,fp);
		fclose(fp);
		free(data);
		}

	return 0;
	}

/* ssn.c */

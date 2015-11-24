/*
 ** txt2spr.c -- DCR 6/11/08
 ** =========
 ** Converts human-readable text to binary springs data file.
 ** =========
 ** Revamped for new springs file format -- SRS 8/4/09
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <rpc/rpc.h>
#include <springs.h>

#define SPR_EXT "spr"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

static void convert(char *infile)
{
	XDR xdrs;
	FILE *fpi,*fpo;
	char outfile[MAXPATHLEN];
	double dTime;
	float fZeroStrainLength,fYoungsModulus,fStressLimit;
	int i,j,nSprings,nPart,iOrder;

	(void) printf("%s\n",infile);
	(void) snprintf(outfile,MAXPATHLEN,"%s.%s",infile,SPR_EXT);
	assert(strcmp(infile,outfile));
	fpi = fopen(infile,"r");
	assert(fpi != NULL);
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	xdrstdio_create(&xdrs,fpo,XDR_ENCODE);
	dTime = 0.0; /* time is not recorded in the text file */
	xdr_double(&xdrs,&dTime);
	(void) fscanf(fpi,"%i",&nPart);
	if (nPart < 1) {
		(void) fprintf(stderr,"%s contains invalid data (nData = %i < 1).\n",infile,nPart);
		exit(1);
		}
	xdr_int(&xdrs,&nPart);
	i = -1; /* pad */
	xdr_int(&xdrs,&i);
	for (i=0;i<nPart;i++) {
		(void) fscanf(fpi,"%i",&nSprings);
		if (nSprings < 0 || nSprings > MAX_NUM_SPRINGS_PER_PARTICLE) {
			(void) fprintf(stderr,"%s contains invalid data (particle %i: nSprings = %i > %i).\n",infile,i,nSprings,MAX_NUM_SPRINGS_PER_PARTICLE);
			exit(1);
		}
	        xdr_int(&xdrs,&nSprings);
		for (j=0;j<nSprings;j++) {
			(void) fscanf(fpi,"%i",&iOrder);
			if (iOrder >= nPart) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: iOrder = %i >= %i).\n",infile,i,j,iOrder,nPart);
				exit(1);
			}
			xdr_int(&xdrs,&iOrder);
			(void) fscanf(fpi,"%f",&fZeroStrainLength);
			if (fZeroStrainLength < 0.0) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fZeroStrainLength = %g < 0).\n",infile,i,j,fZeroStrainLength);
				exit(1);
			}
			xdr_float(&xdrs,&fZeroStrainLength);
			(void) fscanf(fpi,"%f",&fYoungsModulus);
			if (fYoungsModulus < 0.0) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fYoungsModulus = %g < 0).\n",infile,i,j,fYoungsModulus);
				exit(1);
			}
			xdr_float(&xdrs,&fYoungsModulus);
			(void) fscanf(fpi,"%f",&fStressLimit);
			if (fStressLimit < 0.0) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fStressLimit = %g < 0).\n",infile,i,j,fStressLimit);
				exit(1);
			}
			xdr_float(&xdrs,&fStressLimit);
		}
	}
	xdr_destroy(&xdrs);
	(void) fclose(fpo);
	(void) fclose(fpi);
	}

int main(int argc,char *argv[])
{
	int i;

	setbuf(stdout,(char *)NULL);
	if (argc <= 1) {
		(void) fprintf(stderr,"Usage: %s file [ file ... ]\n",argv[0]);
		return 1;
		}
	(void) printf("MAX_NUM_SPRINGS_PER_PARTICLE = %i\n",MAX_NUM_SPRINGS_PER_PARTICLE);
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* txt2spr.c */

/*
** dem2txt.c -- SRS 7/11/11
** =========
** Converts binary DEM data file to human-readable text based upon spr2txt.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <rpc/rpc.h>
/*#include <dem.h>*//*DEBUG*/

#define TXT_EXT "txt"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

static void convert(char *infile)
{
	XDR xdrs;
	FILE *fpi,*fpo;
	char outfile[MAXPATHLEN],cWallsDefined;
	double dTime;
	int i,j,k,nPart,nPE,nWE,iOrder,iWallID;
#if defined(__APPLE__) && defined (__LP64__)
        int liOverlapCounter;
#else
	long int liOverlapCounter;
#endif
	double vShear[3],vnOld[3];

	nWE = 0;
	(void) printf("%s\n",infile);
	(void) snprintf(outfile,MAXPATHLEN,"%s.%s",infile,TXT_EXT);
	assert(strcmp(infile,outfile));
	fpi = fopen(infile,"r");
	assert(fpi != NULL);
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	xdrstdio_create(&xdrs,fpi,XDR_DECODE);
	xdr_double(&xdrs,&dTime);
	if (dTime < 0.0) {
		(void) fprintf(stderr,"%s contains invalid data (dTime = %g < 0.0).\n",infile,dTime);
		exit(1);
		}
	xdr_int(&xdrs,&nPart);
	if (nPart < 1) {
		(void) fprintf(stderr,"%s contains invalid data (nData = %i < 1).\n",infile,nPart);
		exit(1);
		}
	(void) fprintf(fpo,"%i ",nPart);
	xdr_int(&xdrs,&i); /* dummy read */
	xdr_char(&xdrs,&cWallsDefined);
	(void) fprintf(fpo,"%c \n",cWallsDefined);
	if (cWallsDefined != 'P' && cWallsDefined != 'W') {
		(void) fprintf(stderr,"%s contains invalid data (cWallsDefined %c (needs to be P for no !WALLS, W for WALLS)).\n",infile,cWallsDefined);
		exit(1);
		}
	for (i=0;i<nPart;i++) {
		xdr_int(&xdrs,&nPE);
		(void) fprintf(fpo,"%i%s",nPE,
					   nPE == 0 && cWallsDefined == 'P' ? "\n" : " ");
		for (j=0;j<nPE;j++) {
			xdr_int(&xdrs,&iOrder);
			if (iOrder >= nPart) /* just give a warning since pID could be higher than nPart in certain cercumstances (e.g. deleted particles) */
				fprintf(stderr,"WARNING: %s contains suspicious data (particle %i DEM element %i: iOrder = %i >= %i).\n",infile,i,j,iOrder,nPart);
			for (k=0;k<3;k++) xdr_double(&xdrs,&vShear[k]);
			for (k=0;k<3;k++) xdr_double(&xdrs,&vnOld[k]);
			xdr_long(&xdrs,&liOverlapCounter);
			if (liOverlapCounter <= 0) {
				fprintf(stderr,"%s contains invalid data (particle %i DEM element %i: liOverlapCounter = %li <= 0).\n",infile,i,j,(long) liOverlapCounter);
				exit(1);
				}
			fprintf(fpo,"%i %.7e %.7e %.7e %.7e %.7e %.7e %li%s",iOrder,vShear[0],vShear[1],vShear[2],vnOld[0],vnOld[1],vnOld[2],(long) liOverlapCounter,
						   (j == nPE - 1 && cWallsDefined == 'P') ? "\n" : " ");
			}
		if (cWallsDefined == 'W') {
			xdr_int(&xdrs,&nWE);
			(void) fprintf(fpo,"%i%s",nWE,
						   nWE == 0 ? "\n" : " ");
			for (j=0;j<nWE;j++) {
				xdr_int(&xdrs,&iWallID);
				if (iWallID >= nPart) {
					(void) fprintf(stderr,"%s contains invalid data (particle %i DEM wallelement %i: iOrder = %i >= %i).\n",infile,i,j,iWallID,nPart);
					exit(1);
					}
				for (k=0;k<3;k++) xdr_double(&xdrs,&vShear[k]);
				for (k=0;k<3;k++) xdr_double(&xdrs,&vnOld[k]);
				xdr_long(&xdrs,&liOverlapCounter);
				if (liOverlapCounter <= 0) {
					fprintf(stderr,"%s contains invalid data (particle %i DEM wallelement %i: liOverlapCounter = %li <= 0).\n",infile,i,j,(long) liOverlapCounter);
					exit(1);
					}
				fprintf(fpo,"%i %.7e %.7e %.7e %.7e %.7e %.7e %li%s",iWallID,vShear[0],vShear[1],vShear[2],vnOld[0],vnOld[1],vnOld[2],(long) liOverlapCounter,
							   j == nWE - 1 ? "\n" : " ");
				}
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
/*	(void) printf("MAX_NUM_OVERLAPS_PER_PARTICLE = %i\n",MAX_NUM_OVERLAPS_PER_PARTICLE); *//*DEBUG*/
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* dem2txt.c */

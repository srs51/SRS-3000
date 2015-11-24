/*
 ** txt2dem.c -- SRS 8/13/11
 ** =========
 ** Converts human-readable text to binary DEM data file.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <rpc/rpc.h>
/*#include <dem.h>*//*DEBUG*/

#define DEM_EXT "dem"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

static void convert(char *infile)
{
	XDR xdrs;
	FILE *fpi,*fpo;
	char outfile[MAXPATHLEN],cWallsDefined;
	double dTime;
	int i,j,nPart,nPE,nWE,iOrder,iWallID;
#if defined(__APPLE__) && defined (__LP64__)
        int liOverlapCounter;
#else
	long int liOverlapCounter;
#endif
	double vShear[3],vnOld[3];

	nWE = 0;
	(void) printf("%s\n",infile);
	(void) snprintf(outfile,MAXPATHLEN,"%s.%s",infile,DEM_EXT);
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
	(void) fscanf(fpi,"%c",&cWallsDefined); /* skip one space */
	printf("%c\n",cWallsDefined);
	(void) fscanf(fpi,"%c",&cWallsDefined);
	printf("%c\n",cWallsDefined);
	if (cWallsDefined != 'P' && cWallsDefined != 'W') {
		(void) fprintf(stderr,"%s contains invalid data (cWallsDefined %c (needs to be P for no WALLS, W for WALLS)).\n",infile,cWallsDefined);
		exit(1);
		}
	xdr_char(&xdrs,&cWallsDefined);
	for (i=0;i<nPart;i++) {
		(void) fscanf(fpi,"%i",&nPE);
		/* The following check is skipped. We will rely on pkdReadSS() to ensure slots exist for each contact */
		/*
		if (nPE < 0 || nPE > MAX_NUM_OVERLAPS_PER_PARTICLE) {
			(void) fprintf(stderr,"%s contains invalid data (particle %i has %i contacts. Must have 0 <= number of contacts < %i.\n",infile,i,nPE,MAX_NUM_OVERLAPS_PER_PARTICLE);
			exit(1);
		}
		*/
		xdr_int(&xdrs,&nPE);
		for (j=0;j<nPE;j++) {
                        (void) fscanf(fpi,"%i",&iOrder);
                        if (iOrder >= nPart) {
                                (void) fprintf(stderr,"%s contains invalid data (particle %i overlap %i: iOrder = %i >= %i).\n",infile,i,j,iOrder,nPart);
                                exit(1);
                        }
                        xdr_int(&xdrs,&iOrder);

			/* tangential spring read */
                        (void) fscanf(fpi,"%lf",&vShear[0]);  /* error checks for vShear? */
                        xdr_double(&xdrs,&vShear[0]);
                        (void) fscanf(fpi,"%lf",&vShear[1]);
                        xdr_double(&xdrs,&vShear[1]);
                        (void) fscanf(fpi,"%lf",&vShear[2]);
                        xdr_double(&xdrs,&vShear[2]);

			/* previous normal vector read */
                        (void) fscanf(fpi,"%lf",&vnOld[0]);  /* error checks for vnOld? maybe ensure |vnOld| != 0 and warn and normalize if |vnOld| != 1 */
                        xdr_double(&xdrs,&vnOld[0]);
                        (void) fscanf(fpi,"%lf",&vnOld[1]);
                        xdr_double(&xdrs,&vnOld[1]);
                        (void) fscanf(fpi,"%lf",&vnOld[2]);
                        xdr_double(&xdrs,&vnOld[2]);

#if defined(__APPLE__) && defined (__LP64__)
			(void) fscanf(fpi,"%i",&liOverlapCounter);
#else
			(void) fscanf(fpi,"%li",&liOverlapCounter);
#endif
			if (liOverlapCounter <= 0) {
				fprintf(stderr,"%s contains invalid data (particle %i DEM element %i: liOverlapCounter = %li <= 0).\n",infile,i,j,(long) liOverlapCounter);
				exit(1);
				}
			xdr_long(&xdrs,&liOverlapCounter);
			}
		if (cWallsDefined == 'W') {
			(void) fscanf(fpi,"%i",&nWE);
			/* The following check is skipped. We will rely on pkdReadSS() to ensure slots exist for each contact */
			/*
			if (nWE < 0 || nWE > MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i has %i wall contacts. Must have 0 <= number of wall contacts < %i.\n",infile,i,nWE,MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS);
				exit(1);
				}
			*/
			xdr_int(&xdrs,&nWE);
			for (j=0;j<nWE;j++) {
        	                (void) fscanf(fpi,"%i",&iWallID);
                	        if (iOrder >= nPart) {
                        	        (void) fprintf(stderr,"%s contains invalid data (particle %i overlap %i: iWallID = %i >= %i).\n",infile,i,j,iWallID,nPart);
                                	exit(1);
					}
	                        xdr_int(&xdrs,&iWallID);

				/* tangential spring read */
				(void) fscanf(fpi,"%lf",&vShear[0]);  /* error checks for vShear? */
	                        xdr_double(&xdrs,&vShear[0]);
        	                (void) fscanf(fpi,"%lf",&vShear[1]);
                	        xdr_double(&xdrs,&vShear[1]);
                        	(void) fscanf(fpi,"%lf",&vShear[2]);
            		        xdr_double(&xdrs,&vShear[2]);

				/* previous normal vector read */
        	                (void) fscanf(fpi,"%lf",&vnOld[0]);  /* error checks for vnOld? maybe ensure |vnOld| != 0 and warn and normalize if |vnOld| != 1 */
                	        xdr_double(&xdrs,&vnOld[0]);
                        	(void) fscanf(fpi,"%lf",&vnOld[1]);
				xdr_double(&xdrs,&vnOld[1]);
                        	(void) fscanf(fpi,"%lf",&vnOld[2]);
                        	xdr_double(&xdrs,&vnOld[2]);

#if defined(__APPLE__) && defined (__LP64__)
				fscanf(fpi,"%i",&liOverlapCounter);
#else
				fscanf(fpi,"%li",&liOverlapCounter);
#endif
				if (liOverlapCounter <= 0) {
					fprintf(stderr,"%s contains invalid data (particle %i DEM element %i: liOverlapCounter = %li <= 0).\n",infile,i,j,(long) liOverlapCounter);
					exit(1);
					}
				xdr_long(&xdrs,&liOverlapCounter);
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

/* txt2dem.c */

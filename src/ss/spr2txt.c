/*
 ** spr2txt.c -- DCR 6/11/08
 ** =========
 ** Converts binary springs data file to human-readable text.
 ** =========
 ** Revamped for new springs file format -- SRS 8/4/09
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <rpc/rpc.h>
#include <springs.h>

#define TXT_EXT "txt"

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
	(void) fprintf(fpo,"%i\n",nPart);
	xdr_int(&xdrs,&i); /* dummy read */
	for (i=0;i<nPart;i++) {
	        xdr_int(&xdrs,&nSprings);
		(void) fprintf(fpo,"%i%s",nSprings,
			                   nSprings == 0 ? "\n" : " ");
		for (j=0;j<nSprings;j++) {
			xdr_int(&xdrs,&iOrder);
			if (iOrder >= nPart) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: iOrder = %i >= %i).\n",infile,i,j,iOrder,nPart);
				exit(1);
				}
			xdr_float(&xdrs,&fZeroStrainLength);
			if (fZeroStrainLength < 0.) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fZeroStrainLength = %g < 0).\n",infile,i,j,fZeroStrainLength);
				exit(1);
				}
			xdr_float(&xdrs,&fYoungsModulus);
			if (fYoungsModulus < 0.) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fYoungsModulus = %g < 0).\n",infile,i,j,fYoungsModulus);
				exit(1);
				}
			xdr_float(&xdrs,&fStressLimit);
			if (fStressLimit < 0.) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fStressLimit = %g < 0).\n",infile,i,j,fStressLimit);
				exit(1);
				}
			(void) fprintf(fpo,"%i %.7e %.7e %.7e%s",iOrder,fZeroStrainLength,fYoungsModulus,fStressLimit,
						   j == nSprings - 1 ? "\n" : " ");
			/* TEMPORARY READOUT
			printf("%.7e %.7e\n",fYoungsModulus,fStressLimit);
			*/
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

/* spr2txt.c */

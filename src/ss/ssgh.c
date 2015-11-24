/*
 * ssgh.c -- DCR 97-10-21
 * =====
 * Constructs Solar System particle history from genealogy files.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* for isdigit() */
#include <assert.h>
#include <sys/param.h>	/* for MAXPATHLEN */
#include <rpc/rpc.h>	/* for XDR routines */
#include "ss.h"

/*** END OF PREAMBLE ***/

int main(int argc,char *argv[])
{
	extern int getopt(int,char *const *,const char *);
	extern char *optarg;
	extern int optind;

	int get_orig_idx(char *,int,int *);
	int get_curr_idx(char *,int,int *);
	int write_data(char *,int,FILE *);

	FILE *fpout;
	char *mapfile = NULL,outfile[MAXPATHLEN];
	int c,ci=-1,oi=-1;

(void) fprintf(stderr,"THIS CODE IS DEPRECATED -- USE SST INSTEAD\n");
return 1;

	/* Disable stdout buffering */

	setbuf(stdout,(char *) NULL);

	/* Check arguments */

	while ((c = getopt(argc,argv,"m:p:P:")) != EOF)
		switch (c) {
		case 'm':
			mapfile = strdup(optarg);
			break;
		case 'p':
			if (isdigit((int) optarg[0])) ci = atoi(optarg);
			break;
		case 'P':
			if (isdigit((int) optarg[0])) oi = atoi(optarg);
			}

	if ((ci < 0 && oi < 0) || (ci >= 0 && oi >= 0) || (mapfile && ci < 0) ||
		(ci >= 0 && !mapfile) || optind == argc) {
		(void) fprintf(stderr,"Usage: %s ( -m mapfile -p particle# | "
					   "-P particle# ) datfile [ datfile ... ]\n",argv[0]);
		exit(1);
		}

	/*
	 ** Get original index corresponding to requested "current" particle
	 ** index in supplied map file, if applicable.
	 */

	if (mapfile && get_orig_idx(mapfile,ci,&oi))
		exit(1);

	/* Open output file */

	(void) sprintf(outfile,"%i%s",oi,SSH_EXT);

	if (!(fpout = fopen(outfile,"w"))) {
		(void) fprintf(stderr,"Unable to open \"%s\" for writing.\n",outfile);
		exit(1);
		}

	/* Allocate space for new mapfile names */

	mapfile = (char *) realloc(mapfile,MAXPATHLEN);
	assert(mapfile != NULL);

	/* Trace history */

	for (;optind < argc;optind++) {
		(void) printf("%s: ",argv[optind]);
		assert(strlen(argv[optind]) + strlen(MAP_EXT) < MAXPATHLEN);
		(void) sprintf(mapfile,"%s%s",argv[optind],MAP_EXT);
		if (get_curr_idx(mapfile,oi,&ci))
			continue;
		(void) printf("%i --> %i\n",oi,ci);
		(void) write_data(argv[optind],ci,fpout);
		}

	/* All done */

	free((void *) mapfile);

	(void) fclose(fpout);

	return 0;
	}

int get_orig_idx(char *mapfile,int ci,int *oi)
{
	/* try to find "ci" in mapfile, return line number in oi */
	/* error if more than one occurence of ci in mapfile */

	FILE *fp;
	int i,map;

	if (!(fp = fopen(mapfile,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",mapfile);
		return 1;
		}

	*oi = -1;

	for (i = 0;fscanf(fp,"%i",&map) == 1;i++)
		if (map == ci) {
			if (*oi >= 0) {
				(void) fprintf(stderr,"More than one match to index %i found "
							   "in mapfile.\nUse -P option instead.\n",ci);
				return 1;
				}
			*oi = i;
			}

	(void) fclose(fp);

	if (*oi < 0) {
		(void) fprintf(stderr,"Unable to find index %i in mapfile.\n",ci);
		return 1;
		}

	return 0;
	}

int get_curr_idx(char *mapfile,int oi,int *ci)
{
	/* set "ci" to value at line "oi" in "mapfile" */

	FILE *fp;
	int i,map;

	if (!(fp = fopen(mapfile,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",mapfile);
		return 1;
		}

	*ci = -1;

	for (i = 0;fscanf(fp,"%i",&map) == 1;i++)
		if (i == oi) {
			*ci = map;
			break;
			}

	(void) fclose(fp);

	if (*ci == -1) {
		(void) fprintf(stderr,"Unable to reach line %i in mapfile.\n",oi);
		return 1;
		}

	return 0;
	}

int write_data(char *datname,int idx,FILE *fpo)
{
	FILE *fpi;
	XDR xdrs;
	double dum;
	int i,n,pad;

	if (!(fpi = fopen(datname,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",datname);
		return 1;
		}

	xdrstdio_create(&xdrs,fpi,XDR_DECODE);

	(void) xdr_double(&xdrs,&dum); /* time */
	(void) xdr_int(&xdrs,&n);
	(void) xdr_int(&xdrs,&pad);

	if (idx >= n) {
		(void) fprintf(stderr,"Unable to find particle %i.\n",idx);
		xdr_destroy(&xdrs);
		(void) fclose(fpi);
		return 1;
		}

	if (fseek(fpi,idx*sizeof(SSDATA),SEEK_CUR)) {
		(void) fprintf(stderr,"Unable to seek to particle %i data.\n",idx);
		xdr_destroy(&xdrs);
		(void) fclose(fpi);
		return 1;
		}

	(void) fprintf(fpo,"%.16e",dum);

	for (i=0;i<11;i++) {
		(void) xdr_double(&xdrs,&dum);
		(void) fprintf(fpo," %.16e",dum);
		}

	(void) fprintf(fpo,"\n");

	xdr_destroy(&xdrs);
	(void) fclose(fpi);
	return 0;
	}

/* ssgh.c */

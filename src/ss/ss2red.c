/*
 ** ss2red.c -- DCR 6/17/08
 ** =======-
 ** Converts Solar System binary data to reduced form.
 */

#include <ssio.h>
#include <vector.h>

#define RED_EXT ".r"

static void convert(char *infile)
{
	SSIO ssio_i,ssio_o;
	SSHEAD head;
	SSDATA d;
	SSRDATA dr;
	char outfile[MAXPATHLEN];
	int i;

	(void) printf("%s: ",infile);
	(void) sprintf(outfile,"%s%s",infile,RED_EXT); /*DEBUG check for overflow!*/
	if (ssioOpen(infile,&ssio_i,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\" for reading.\n",infile);
		return;
		}
	if (ssioOpen(outfile,&ssio_o,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\" for writing.\n",outfile);
		(void) ssioClose(&ssio_i);
		return;
		}
	if (ssioHead(&ssio_i,&head)) {
		(void) fprintf(stderr,"Corrupt header.\n");
		goto finish;
		}
	(void) printf("time = %e, n_data = %i\n",head.time,head.n_data);
	if (head.n_data <= 0) {
		(void) fprintf(stderr,"Invalid input data format.\n");
		goto finish;
		}
	switch(head.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"File already in reduced format.\n");
		goto finish;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",head.iMagicNumber);
		goto finish;
		}
	head.iMagicNumber = SSIO_MAGIC_REDUCED;
	if (ssioHead(&ssio_o,&head)) {
		(void) fprintf(stderr,"Error writing header.\n");
		goto finish;
		}
	for (i=0;i<head.n_data;i++) {
		if (ssioData(&ssio_i,&d)) {
			(void) fprintf(stderr,"Corrupt data.\n");
			goto finish;
			}
		ssioStandardToReduced(&d,&dr);
		if (ssioDataReduced(&ssio_o,&dr)) {
			(void) fprintf(stderr,"Error writing data.\n");
			goto finish;
			}
		}
 finish:
	(void) ssioClose(&ssio_o);
	(void) ssioClose(&ssio_i);
	}

int main(int argc,char *argv[])
{
	int i;

	setbuf(stdout,(char *)NULL);
	if (argc <= 1) {
		(void) fprintf(stderr,"Usage: %s file [ file ... ]\n",argv[0]);
		exit(1);
		}
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* ss2red.c */

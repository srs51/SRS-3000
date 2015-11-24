/*
 ** random.c -- DCR 3/9/12
 ** ========
 **
 ** Implementation of NR(3e) uniform & Gaussian RNGs.
 **
 ** Compare random.c in the pkdgrav source code.
 */

#include <stdio.h>
#include <math.h>
#include <random.h>

static Ullong ullState = 4101842887655102017LL; /* used by randSeedGenerator() -- don't use this value as the seed! */

static Ullong rand_uniform_int64(void) /* Ranq1.int64() in NR(3e) */
/*
** Recommended generator for everyday use.  The period is ~1.8e19.
*/
{
	ullState ^= ullState >> 21;
	ullState ^= ullState << 35;
	ullState ^= ullState >> 4;

	return ullState*2685821657736338717LL;
	}

double randUniform(void) /* Ranq1.doub() in NR(3e) */
{
	return 5.42101086242752217e-20*rand_uniform_int64();
	}

double randGaussian(void) /* Normaldev.dev() in NR(3e) */
{
	double u,v,x,y,q;

	do {
		u = randUniform();
		v = 1.7156*(randUniform() - 0.5);
		x = u - 0.449871;
		y = fabs(v) + 0.386595;
		q = x*x + y*(0.19600*y - 0.25472*x);
		} while (q > 0.27597 && (q > 0.27846 || v*v > -4.0*log(u)*u*u));

	return v/u; /* zero mean, unit standard deviation */
	}

void randSeedGenerator(int iSeed) /* based on Ranq1 constructor in NR(3e) */
{
	Ullong ullSeed = (Ullong) iSeed;

	ullState ^= ullSeed;
	ullState = rand_uniform_int64();
	}

/* for testing only */

/*#define TEST*/
#ifdef TEST

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

int main(int argc,char *argv[])
{
	int iSeed = (int) time(NULL) % getpid() + getppid();

	randSeedGenerator(iSeed);
	printf("Uniform random deviate = %g\n",randUniform());
	printf("Gaussian deviate, unit std dev = %g\n",randGaussian());

	return 0;
	}

#endif /* TEST */

/* random.c */

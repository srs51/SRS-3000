/*
 ** NULL mdl.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mdl.h"


#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2


#define MDL_DEFAULT_BYTES		4096
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

void _srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut)
{
	return;
	}


double mdlCpuTimer(MDL mdl)
{
#ifndef _CRAYMPP
	struct rusage ru;

	getrusage(0,&ru);
	return((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
#else
	return 0.0;
#endif
	}

/* 
 * MDL debug and Timer functions 
 */
#define MDLPRINTF_STRING_MAXLEN 256
void mdlprintf( MDL mdl, const char *format, ... )
{
     static char ach[MDLPRINTF_STRING_MAXLEN];
     va_list args;

     if (mdl->bDiag) {	
         va_start( args, format);
         vsnprintf( ach, MDLPRINTF_STRING_MAXLEN, format, args);
         mdlDiag( mdl, ach);
         va_end( args);
         }
}

#ifdef MDLDEBUG
void mdldebug( MDL mdl, const char *format, ... )
{
     static char ach[MDLPRINTF_STRING_MAXLEN];
     va_list args;

     if (mdl->bDiag) {	
         va_start( args, format);
	 vsnprintf( ach, MDLPRINTF_STRING_MAXLEN, format, args);
	 mdlDiag( mdl, ach);
	 va_end( args);
         }
}
#endif

#ifdef MDLTIMER
void mdlZeroTimer(MDL mdl, mdlTimer *t)
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  t->wallclock = tv.tv_sec + 1e-6*(double) tv.tv_usec;

  struct rusage ru;
  getrusage(0,&ru);
  t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec;
  t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
}

void mdlGetTimer(MDL mdl, mdlTimer *t0, mdlTimer *t)
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  t->wallclock = tv.tv_sec + 1e-6*(double) tv.tv_usec - t0->wallclock;

  struct rusage ru;
  getrusage(0,&ru);
  t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec - t0->cpu;
  t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec - t0->system;
}

void mdlPrintTimer(MDL mdl,char *message, mdlTimer *t0) 
{
  mdlTimer lt;

  if (mdl->bDiag) {	
      mdlGetTimer(mdl,t0,&lt);
      mdlprintf(mdl,"%s %f %f %f\n",message,lt.wallclock,lt.cpu,lt.system);
      }
}
#endif

int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{
	MDL mdl;
	int i,nThreads,bDiag,bThreads;
	char *p,ach[256],achDiag[256];

	*pmdl = NULL;
	mdl = malloc(sizeof(struct mdlContext));
	assert(mdl != NULL);
	/*
	 ** Set default "maximums" for structures. These are NOT hard
	 ** maximums, as the structures will be realloc'd when these
	 ** values are exceeded.
	 */
	mdl->nMaxServices = MDL_DEFAULT_SERVICES;
	mdl->nMaxInBytes = MDL_DEFAULT_BYTES;
	mdl->nMaxOutBytes = MDL_DEFAULT_BYTES;
	mdl->nMaxCacheIds = MDL_DEFAULT_CACHEIDS;
	/*
	 ** Now allocate the initial service slots.
	 */
	mdl->psrv = malloc(mdl->nMaxServices*sizeof(SERVICE));
	assert(mdl->psrv != NULL);
	/*
	 ** Initialize the new service slots.
	 */
	for (i=0;i<mdl->nMaxServices;++i) {
		mdl->psrv[i].p1 = NULL;
		mdl->psrv[i].nInBytes = 0;
		mdl->psrv[i].nOutBytes = 0;
		mdl->psrv[i].fcnService = NULL;
		}
	/*
	 ** Provide a 'null' service for sid = 0, so that stopping the 
	 ** service handler is well defined!
	 */
	mdl->psrv[0].p1 = NULL;
	mdl->psrv[0].nInBytes = 0;
	mdl->psrv[0].nOutBytes = 0;
	mdl->psrv[0].fcnService = _srvNull;
	/*
	 ** Allocate service buffers.
	 */
	mdl->pszIn = malloc(mdl->nMaxInBytes);
	assert(mdl->pszIn != NULL);
	mdl->pszOut = malloc(mdl->nMaxOutBytes);
	assert(mdl->pszOut != NULL);
	/*
	 ** Allocate initial cache spaces.
	 */
	mdl->cache = malloc(mdl->nMaxCacheIds*sizeof(CACHE));
	assert(mdl->cache != NULL);
	/*
	 ** Initialize caching spaces.
	 */
	for (i=0;i<mdl->nMaxCacheIds;++i) {
		mdl->cache[i].iType = MDL_NOCACHE;
		}
	/*
	 ** Do some low level argument parsing for number of threads, and
	 ** diagnostic flag!
	 */
	bDiag = 0;
	bThreads = 0;
	i = 1;
	while (argv[i]) {
		if (!strcmp(argv[i],"-sz") && !bThreads) {
			++i;
			if (argv[i]) {
				nThreads = atoi(argv[i]);
				bThreads = 1;
				}
			}
		if (!strcmp(argv[i],"+d") && !bDiag) {
			p = getenv("MDL_DIAGNOSTIC");
			if (!p) p = getenv("HOME");
			if (!p) sprintf(ach,"/tmp");
			else sprintf(ach,"%s",p);
			bDiag = 1;
			}
		++i;
		}
	nThreads = 1;
	mdl->bDiag = bDiag;
	mdl->nThreads = nThreads;
	*pmdl = mdl;
	/*
	 ** A unik!
	 */
	mdl->idSelf = 0;
	if (mdl->bDiag) {
		char *tmp = strrchr(argv[0],'/');
		if (!tmp) tmp = argv[0];
		else ++tmp;
		sprintf(achDiag,"%s/%s.%d",ach,tmp,mdl->idSelf);
		mdl->fpDiag = fopen(achDiag,"w");
		assert(mdl->fpDiag != NULL);
		}
	return(nThreads);
	}


void mdlFinish(MDL mdl)
{
	/*
	 ** Close Diagnostic file.
	 */
	if (mdl->bDiag) {
		fclose(mdl->fpDiag);
		}
	/*
	 ** Deregister from PVM and deallocate storage.
	 */
	free(mdl->psrv);
	free(mdl->pszIn);
	free(mdl->pszOut);
	free(mdl->cache);
	free(mdl->atid);
	free(mdl);
	}


/*
 ** This function returns the number of threads in the set of 
 ** threads.
 */
int mdlThreads(MDL mdl)
{
	return(mdl->nThreads);
	}


/*
 ** This function returns this threads 'id' number within the specified
 ** MDL Context. Parent thread always has 'id' of 0, where as children
 ** have 'id's ranging from 1..(nThreads - 1).
 */
int mdlSelf(MDL mdl)
{
	return(mdl->idSelf);
	}


/*
 ** This is a tricky function. It initiates a bilateral transfer between
 ** two threads. Both threads MUST be expecting this transfer. The transfer
 ** occurs between idSelf <---> 'id' or 'id' <---> idSelf as seen from the
 ** opposing thread. It is designed as a high performance non-local memory
 ** swapping primitive and implementation will vary in non-trivial ways 
 ** between differing architectures and parallel paradigms (eg. message
 ** passing and shared address space). A buffer is specified by 'pszBuf'
 ** which is 'nBufBytes' in size. Of this buffer the LAST 'nOutBytes' are
 ** transfered to the opponent, in turn, the opponent thread transfers his
 ** nBufBytes to this thread's buffer starting at 'pszBuf'.
 ** If the transfer completes with no problems the function returns 1.
 ** If the function returns 0 then one of the players has not received all
 ** of the others memory, however he will have successfully transfered all
 ** of his memory.
 */
int mdlSwap(MDL mdl,int id,size_t nBufBytes,void *vBuf,size_t nOutBytes,
			size_t *pnSndBytes,size_t *pnRcvBytes)
{
	assert(0);
	return 0;
	}


void mdlDiag(MDL mdl,char *psz)
{
	if (mdl->bDiag) {	
		fputs(psz,mdl->fpDiag);
		fflush(mdl->fpDiag);
		}
	}


void mdlAddService(MDL mdl,int sid,void *p1,
				   void (*fcnService)(void *,void *,int,void *,int *),
				   int nInBytes,int nOutBytes)
{
	int i,nMaxServices;

	assert(sid > 0);
	if (sid >= mdl->nMaxServices) {
		/*
		 ** reallocate service buffer, adding space for 8 new services
		 ** including the one just defined.
		 */
		nMaxServices = sid + 9;
		mdl->psrv = realloc(mdl->psrv,nMaxServices*sizeof(SERVICE));
		assert(mdl->psrv != NULL);
		/*
		 ** Initialize the new service slots.
		 */
		for (i=mdl->nMaxServices;i<nMaxServices;++i) {
			mdl->psrv[i].p1 = NULL;
			mdl->psrv[i].nInBytes = 0;
			mdl->psrv[i].nOutBytes = 0;
			mdl->psrv[i].fcnService = NULL;
			}
		mdl->nMaxServices = nMaxServices;
		}
	/*
	 ** Make sure the service buffers are big enough!
	 */
	if (nInBytes > mdl->nMaxInBytes) {
		mdl->pszIn = realloc(mdl->pszIn,nInBytes);
		assert(mdl->pszIn != NULL);
		mdl->nMaxInBytes = nInBytes;
		}
	if (nOutBytes > mdl->nMaxOutBytes) {
		mdl->pszOut = realloc(mdl->pszOut,nOutBytes);
		assert(mdl->pszOut != NULL);
		mdl->nMaxOutBytes = nOutBytes;
		}
	mdl->psrv[sid].p1 = p1;
	mdl->psrv[sid].nInBytes = nInBytes;
	mdl->psrv[sid].nOutBytes = nOutBytes;
	mdl->psrv[sid].fcnService = fcnService;
	}


void mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes)
{
	assert(0);
	}


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes)
{
	assert(0);
	}


void mdlHandler(MDL mdl)
{
	assert(0);
	}


#define BILLION				1000000000

/*
 ** Special MDL memory allocation functions for allocating memory 
 ** which must be visible to other processors thru the MDL cache 
 ** functions.
 ** mdlMalloc() is defined to return a pointer to AT LEAST iSize bytes 
 ** of memory. This pointer will be passed to either mdlROcache or 
 ** mdlCOcache as the pData parameter.
 ** For PVM and most machines these functions are trivial, but on the 
 ** T3D and perhaps some future machines these functions are required.
 */
void *mdlMalloc(MDL mdl,int iSize)
{	
	return(malloc(iSize));
	}


void mdlFree(MDL mdl,void *p)
{
	free(p);
	}


/*
 ** Common initialization for all types of caches.
 */
CACHE *CacheInitialize(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int i,nMaxCacheIds;

	/*
	 ** Allocate more cache spaces if required!
	 */
	assert(cid >= 0);
	if (cid >= mdl->nMaxCacheIds) {
		/*
		 ** reallocate cache spaces, adding space for 2 new cache spaces
		 ** including the one just defined.
		 */
		nMaxCacheIds = cid + 3;
		mdl->cache = realloc(mdl->cache,nMaxCacheIds*sizeof(CACHE));
		assert(mdl->cache != NULL);
		/*
		 ** Initialize the new cache slots.
		 */
		for (i=mdl->nMaxCacheIds;i<nMaxCacheIds;++i) {
			mdl->cache[i].iType = MDL_NOCACHE;
			}
		mdl->nMaxCacheIds = nMaxCacheIds;
		}
	c = &mdl->cache[cid];
	assert(c->iType == MDL_NOCACHE);
	c->pData = pData;
	c->iDataSize = iDataSize;
	c->nData = nData;
	c->nAccess = 0;
	c->nAccHigh = 0;
	c->nMiss = 0;
	c->nColl = 0;
	c->nMin = 0;
	return(c);
	}


/*
 ** Initialize a Read-Only caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_ROCACHE;
	c->init = NULL;
	c->combine = NULL;
	}


/*
 ** Initialize a Combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
{
	CACHE *c;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_COCACHE;
	c->init = init;
	c->combine = combine;
	}


void mdlFinishCache(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];

	/*
	 ** Free up storage and finish.
	 */
	c->iType = MDL_NOCACHE;
	}

void mdlCacheCheck(MDL mdl)
{
}

void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];

	++c->nAccess;
	assert(id == mdl->idSelf);
	return(&c->pData[iIndex*c->iDataSize]);
	}


void mdlRelease(MDL mdl,int cid,void *p)
{
	}


double mdlNumAccess(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];

	return(c->nAccHigh*1e9 + c->nAccess);
	}


double mdlMissRatio(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	double dAccess = c->nAccHigh*1e9 + c->nAccess;
	
	if (dAccess > 0.0) return(c->nMiss/dAccess);
	else return(0.0);
	}


double mdlCollRatio(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	double dAccess = c->nAccHigh*1e9 + c->nAccess;

	if (dAccess > 0.0) return(c->nColl/dAccess);
	else return(0.0);
	}


double mdlMinRatio(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	double dAccess = c->nAccHigh*1e9 + c->nAccess;

	if (dAccess > 0.0) return(c->nMin/dAccess);
	else return(0.0);
	}

#ifndef _TIMER_INCLUDED_
#define _TIMER_INCLUDED_

#ifdef WIN32
// TODO
#else
#include <time.h>
typedef struct timespec TIMEVAL;
#endif

#if defined(__cplusplus)
extern "C"
{
#endif

void getTime(TIMEVAL *t);
long getDiffNanosecs(TIMEVAL *first, TIMEVAL *second);

#if defined(__cplusplus)
}
#endif

#endif

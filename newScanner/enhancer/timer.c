//copyright (c) in silico Labs, LLC 2008

/*
 * This is a timer interface for a nanosecond resolution
 * timer that works on both windows and unix.
 *
 * use as in:
 *
 *   TIMEVAL time1, time2;
 *   getTime(&time1);
 *   < do something that takes awhile >
 *   getTime(&time2);
 *   cout << "That took: " << setprecision(3) <<
 *           (getDiffNanosecs(&time1, &time2)/1000000000.0) <<
 *           " sec." << endl;
 *
 */

#include "timer.h"

void getTime(TIMEVAL *t) {
#ifdef WIN32
		// TODO
#else
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, t);
#endif
}

long getDiffNanosecs(TIMEVAL *first, TIMEVAL *second)
{
#ifdef WIN32
		return 0;
#else
		return (second->tv_sec - first->tv_sec) * 1000000000 +
				(second->tv_nsec - first->tv_nsec);
#endif
}

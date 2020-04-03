#ifndef TIMERS_H
#define TIMERS_H

#include <time.h>
#include <sys/time.h>


/**
 * @brief Wrapper on timeval struct
 */
typedef struct timer{
    struct timeval start;
    struct timeval end;
}Timer;

/**
 * @brief monitoring function (init start)
 */
void tic(Timer* inTimer){
    gettimeofday(&inTimer->start, NULL);
}

/**
 * @brief monitoring function (init end)
 */
void tac(Timer* inTimer){
    gettimeofday(&inTimer->end, NULL);
}

/**
 * @brief monitoring function (return elapsed time in micro second)
 */
long int get_elapsed(Timer* inTimer){
    return ((inTimer->end.tv_sec * 1000000 + inTimer->end.tv_usec)
            - (inTimer->start.tv_sec * 1000000 + inTimer->start.tv_usec));
}

/**
 * @brief monitoring function (print elapsed)
 */
void print_elapsed(Timer* inTimer){
    long int elapsed = get_elapsed(inTimer);
    printf("Elapsed : %ld us (or %f seconds)\n",elapsed,elapsed/1000000.0);
}

/**
 * @brief monitoring function (print difference :: First-Second)
 */
void print_difference_elapsed(Timer* inTimer1,Timer*inTimer2){
    long int diff = get_elapsed(inTimer1)-get_elapsed(inTimer2);
    printf("Timer Difference : %ld us (or %f seconds)\n",diff,(double)diff/1000000.0);
}

#endif

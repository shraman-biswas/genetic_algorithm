#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* upper and lower ascii character limts */
#define UPPER 126
#define LOWER 32

/* genome structure */
typedef struct __genome_t {
	unsigned int err;
	char *wts;
} genome_t;

/* population structure */
typedef struct __pop_t{
	bool found;
	int size, gene_cnt, tourn_size;
	float c_rate, m_rate;
	char *pwt1, *pwt2, *cwt1, *cwt2;
	const char *target;
	genome_t *pool, *best;
} pop_t;

/* sign table for mutation operation */
static const int sign[] = {-1, 1};

pop_t *ga_create(
	const char *target,
	const int gene_cnt,
	const int size,
	const double c_rate,
 	const double m_rate,
	const int tourn_size);
void ga_destroy(pop_t *const pop);
void ga_epoch(pop_t *const pop);
genome_t *ga_run(pop_t *const pop, const int num_gen);
void disp_genome(const pop_t *const pop, const genome_t *const g);

#endif
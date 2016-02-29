#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define UPPER 90
#define LOWER 65

typedef struct __genome_t {
	unsigned int err;
	char *wts;
} genome_t;

typedef struct __pop_t{
	bool done;
	int size, gene_cnt, tourn_size;
	float c_rate, m_rate;
	char *pwt1, *pwt2, *cwt1, *cwt2;
	const char *target;
	genome_t *pool;
} pop_t;

pop_t *ga_create(
	const char *target,
	const int gene_cnt,
	const int size,
	const double c_rate,
 	const double m_rate,
	const int tourn_size);
void ga_destroy(pop_t *const pop);
void ga_epoch(pop_t *const pop);
bool ga_run(pop_t *const pop, const int num_gen);

#endif
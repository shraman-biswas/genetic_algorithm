#include "genetic_algorithm.h"

/*----------------------------------------------------------------------------*/
/* internal genetic algorithm functions                                       */
/*----------------------------------------------------------------------------*/

/* compare 2 genomes and return fitter genome with lower error */
static int cmp_genome(const void *a, const void *b)
{
	const genome_t *g1=NULL, *g2=NULL;
	g1 = (const genome_t *)a;
	g2 = (const genome_t *)b;
	return g1->err - g2->err;
}

/* store best genome with least error */
static void best_genome(const pop_t *const pop, const genome_t *const g)
{
	int i;
	for (i=0; i < pop->gene_cnt; ++i)
		pop->best->wts[i] = g->wts[i];
	pop->best->err = g->err;
}

/* calculate error of entire population and determine fittest genome */
static void pop_error(pop_t *const pop)
{
	int i, j, err_gene, err_min=INT_MAX;
	genome_t *g=NULL;
	for (i=0; i < pop->size; ++i) {
		/* clear initial genome error */
		pop->pool[i].err = 0;
		/* calculate genome error */
		for (j=0; j < pop->gene_cnt; ++j) {
			err_gene = pop->pool[i].wts[j] - pop->target[j];
			pop->pool[i].err += err_gene * err_gene;
		}
		/* find fittest genome with least error */
		if (pop->pool[i].err <= err_min) {
			g = &pop->pool[i];
			/* break if no error */
			if (pop->pool[i].err == 0) {
				/* save final fittest genome */
				best_genome(pop, g);
				pop->found = true;
				return;
			}
			err_min = pop->pool[i].err;
		}
	}
	/* update fittest genome after entire population search */
	best_genome(pop, g);
}

/* select 2 fit parents from population pool */
static void tourn_select(
	const pop_t *const pop,
	char *const pwt1,
	char *const pwt2)
{
	int i;
	genome_t *p1, *p2, *tmp;
	/* select 2 random parents from population */
	p1 = &pop->pool[rand() % pop->size];
	p2 = &pop->pool[rand() % pop->size];
	/* pick the fitter parent as parent 1 */
	if (p1->err > p2->err) {
		tmp = p1;
		p1 = p2;
		p2 = tmp;
	}
	for (i=0; i < pop->tourn_size-2; ++i) {
		/* pick another random parent from population */
		tmp = &pop->pool[rand() % pop->size];
		/* set as parent 1 if new parent is fitter else as parent 2 */
		if (tmp->err < p1->err) {
			p1 = tmp;
		} else if (tmp->err < p2->err) {
			p2 = tmp;
		}
	}
	/* copy genes of parents */
	for (i=0; i < pop->gene_cnt; ++i) {
		pwt1[i] = p1->wts[i];
		pwt2[i] = p2->wts[i];
	}
}

/* crossover 2 parents to produce 2 children */
static void crossover(
	const pop_t *const pop,
	const char *const pwt1,
	const char *const pwt2,
	char *const cwt1,
	char *const cwt2)
{
	int i, cross_pt;
	/* return parents if crossover rate probabilty does not suffice */
	if ((float)rand() / (float)RAND_MAX > pop->c_rate) {
		for (i=0; i < pop->gene_cnt; ++i) {
			cwt1[i] = pwt1[i];
			cwt2[i] = pwt2[i];
		}
		return;
	}
	/* calculate crossover point */
	cross_pt = rand() % pop->gene_cnt;
	/* copy genes upto crossover point from parents to children */
	for (i=0; i<cross_pt; ++i) {
		cwt1[i] = pwt1[i];
		cwt2[i] = pwt2[i];
	}
	/* swap genes from crossover point from parents to children */
	for (i=cross_pt; i < pop->gene_cnt; ++i) {
		cwt1[i] = pwt2[i];
		cwt2[i] = pwt1[i];
	}
}

/* mutate child to produce mutated child */
static void mutate(const pop_t *const pop, char *const wt)
{
	int i;
	/* loop over every gene of genome */
	for (i=0; i < pop->gene_cnt; ++i) {
		/* skip gene if mutation rate probability does not suffice */
		if ((float)rand() / (float)RAND_MAX > pop->m_rate)
			continue;
		/* check gene mutation bounds and apply mutation */
		if ((LOWER < wt[i]) && (wt[i] < UPPER)) {
			wt[i] += sign[rand() % 2];
		} else if (wt[i] <= LOWER) {
			wt[i] += 1;
		} else {
			wt[i] -= 1;
		}
	}
}

/* allocate memory for genetic algorithm */
static void init(pop_t *const pop)
{
	/* allocate memory for parents and children genes */
	pop->pwt1 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->pwt2 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->cwt1 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->cwt2 = (char *)calloc(pop->gene_cnt, sizeof(char));
	/* allocate memory for population pool */
	pop->pool = (genome_t *)calloc(pop->size, sizeof(genome_t));
	/* allocate memory for best genome */
	pop->best = (genome_t *)malloc(sizeof(genome_t));
	/* allocate memory for best genome genes */
	pop->best->wts = (char *)calloc(pop->gene_cnt, sizeof(char));
}

/*----------------------------------------------------------------------------*/
/* external genetic algorithm functions                                       */
/*----------------------------------------------------------------------------*/

/* create genetic algorithm */
pop_t *ga_create(
	const char *target,
	const unsigned int gene_cnt,
	const unsigned int size,
	const double c_rate,
	const double m_rate,
	const unsigned int tourn_size)
{
	int i, j;
	/* create population */
	pop_t *pop = (pop_t *)malloc(sizeof(pop_t));
	if (pop == NULL) {
		perror("population could not be created!");
		exit(EXIT_FAILURE);
	}
	/* set best genome found state to false */
	pop->found = false;
	/* get parameters */
	pop->target = target;			/* target string */
	pop->gene_cnt = gene_cnt;		/* number of genes */
	pop->size = (size % 2) ? size-1 : size;	/* population size */
	pop->c_rate = c_rate;			/* crossover rate */
	pop->m_rate = m_rate;			/* mutation rate */
	pop->tourn_size = tourn_size;		/* tournament size */
	/* allocate memory */
	init(pop);
	/* initialize random number generator */
	srand(time(NULL));
	/* create entire population pool */
	for (i=0; i<size; ++i) {
		/* create genes of each genome */
		pop->pool[i].wts = (char *)calloc(gene_cnt, sizeof(char));
		for (j=0; j<gene_cnt; ++j)
			/* initialize each genome with random gene values */
			pop->pool[i].wts[j] = rand() %
				(UPPER - LOWER + 1) + LOWER;
	}
	return pop;
}

/* destroy genetic algorithm and deallocate memory */
void ga_destroy(pop_t *const pop)
{
	int i;
	for (i=0; i< pop->size; ++i)
		free(pop->pool[i].wts);		/* free each genome */
	free(pop->pool);			/* free population pool */
	free(pop->pwt1);			/* free parent 1 genes */
	free(pop->pwt2);			/* free parent 2 genes */
	free(pop->cwt1);			/* free child 1 genes */
	free(pop->cwt2);			/* free child 2 genes */
	free(pop);
}

/* one population life cycle */
void ga_epoch(pop_t *const pop)
{
	char *wts[pop->size];
	int i, j;
	/* create new temporary population pool */
	for (i=0; i < pop->size; ++i)
		wts[i] = (char *)calloc(pop->size, sizeof(char));
	/* calculate error of current population and determine fittest genome */
	pop_error(pop);
	/* sort current population pool in ascending order of error */
	qsort(pop->pool, pop->size, sizeof(genome_t), cmp_genome);
	/* loop over current population pool in steps of 2 */
	for (i=0; i < pop->size; i+=2) {
		/* select 2 fit parents by tournament selection */
		tourn_select(pop, pop->pwt1, pop->pwt2);
		/* crossover both parents to produce 2 children */
		crossover(pop, pop->pwt1, pop->pwt2, pop->cwt1, pop->cwt2);
		/* mutate both children */
		mutate(pop, pop->cwt1);
		mutate(pop, pop->cwt2);
		/* copy children genomes into temporary population pool */
		for (j=0; j < pop->gene_cnt; ++j) {
			wts[i][j] = pop->cwt1[j];
			wts[i+1][j] = pop->cwt2[j];
		}
	}
	/* copy temporary population pool into current population pool */
	for (i=0; i < pop->size; ++i) {
		for (j=0; j < pop->gene_cnt; ++j)
			pop->pool[i].wts[j] = wts[i][j];
		/* free temporary population pool genomes */
		free(wts[i]);
	}
}

/* run genetic algorithm for several generations */
genome_t *ga_run(pop_t *const pop, const unsigned int num_gen)
{
	int i;
	/* loop over all generations */
	for (i=0; i<num_gen; ++i) {
		/* run genetic algorithm for 1 generation */
		ga_epoch(pop);
		/* display current best genome and error */
		printf("\rbest:   [ %s ] error: %d generation: %d",
			pop->best->wts, pop->best->err, i+1);
		/* terminate loop if best genome with zero error is found */
		if (pop->found == true) {
			/* return best genome */
			return pop->best;
		}
	}
	/* return NULL if best genome with zero error was not found */
	return NULL;
}

/*---------------------------------------------------------------------------*/
/* external genetic algorithm helper functions                               */
/*---------------------------------------------------------------------------*/

/* display genome string and error value */
inline void disp_genome(const pop_t *const pop, const genome_t *const g)
{
	printf("genome: [ %s ] error: %d\n", g->wts, g->err);
}

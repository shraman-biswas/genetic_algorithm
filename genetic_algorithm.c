#include "genetic_algorithm.h"

/*---------------------------------------------------------------------------*/
/* internal genetic algorithm functions                                      */
/*---------------------------------------------------------------------------*/

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
		pop->pool[i].err = 0;
		for (j=0; j < pop->gene_cnt; ++j) {
			err_gene = pop->pool[i].wts[j] - pop->target[j];
			pop->pool[i].err += err_gene * err_gene;
		}
		if (pop->pool[i].err <= err_min) {
			g = &pop->pool[i];
			if (pop->pool[i].err == 0) {
				best_genome(pop, g);
				pop->found = true;
				return;
			}
			err_min = pop->pool[i].err;
		}
	}
	best_genome(pop, g);
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
	if ((float)rand() / (float)RAND_MAX > pop->c_rate) {
		for (i=0; i < pop->gene_cnt; ++i) {
			cwt1[i] = pwt1[i];
			cwt2[i] = pwt2[i];
		}
		return;
	}
	cross_pt = rand() % pop->gene_cnt;
	for (i=0; i<cross_pt; ++i) {
		cwt1[i] = pwt1[i];
		cwt2[i] = pwt2[i];
	}
	for (i=cross_pt; i< pop->gene_cnt; ++i) {
		cwt1[i] = pwt2[i];
		cwt2[i] = pwt1[i];
	}
}

/* mutate child to produce mutated child */
static void mutate(const pop_t *const pop, char *const wt)
{
	int i;
	for (i=0; i < pop->gene_cnt; ++i) {
		if ((float)rand() / (float)RAND_MAX > pop->m_rate)
			continue;
		if ((LOWER < wt[i]) && (wt[i] < UPPER)) {
			wt[i] += sign[rand() % 2];
		} else if (wt[i] <= LOWER) {
			wt[i] += 1;
		} else {
			wt[i] -= 1;
		}
	}
}

/* select 2 fit parents from population pool */
static void tourn_select(
	const pop_t *const pop,
	char *const pwt1,
	char *const pwt2)	
{
	int i;
	genome_t *p1, *p2, *tmp;
	p1 = &pop->pool[rand() % pop->size];
	p2 = &pop->pool[rand() % pop->size];
	if (p1->err > p2->err) {
		tmp = p1;
		p1 = p2;
		p2 = tmp;
	}
	for (i=0; i < pop->tourn_size-2; ++i) {
		tmp = &pop->pool[rand() % pop->size];
		if (tmp->err < p1->err) {
			p1 = tmp;
		} else if (tmp->err < p2->err) {
			p2 = tmp;
		}
	}
	for (i=0; i < pop->gene_cnt; ++i) {
		pwt1[i] = p1->wts[i];
		pwt2[i] = p2->wts[i];
	}
}

/* allocate memory for genetic algorithm */
static void init(pop_t *const pop)
{
	pop->pwt1 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->pwt2 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->cwt1 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->cwt2 = (char *)calloc(pop->gene_cnt, sizeof(char));
	pop->pool = (genome_t *)calloc(pop->size, sizeof(genome_t));
	pop->best = (genome_t *)malloc(sizeof(genome_t));
	pop->best->wts = (char *)calloc(pop->gene_cnt, sizeof(char));	
}

/*---------------------------------------------------------------------------*/
/* external genetic algorithm functions                                      */
/*---------------------------------------------------------------------------*/

/* create genetic algorithm */
pop_t *ga_create(
	const char *target,
	const int gene_cnt,
	const int size,
	const double c_rate,
	const double m_rate,
	const int tourn_size)
{
	int i, j;
	pop_t *pop = (pop_t *)malloc(sizeof(pop_t));
	if (pop == NULL) {
		perror("population could not be created!");
		exit(EXIT_FAILURE);
	}
	pop->found = false;
	/* get parameters */
	pop->target = target;						/* target string */
	pop->gene_cnt = gene_cnt;					/* number of genes */
	pop->size = (size % 2) ? size-1 : size;		/* population size */
	pop->c_rate = c_rate;						/* crossover rate */
	pop->m_rate = m_rate;						/* mutation rate */
	pop->tourn_size = tourn_size;				/* tournament size */
	/* allocate memory */
	init(pop);
	/* initialize random number generator */
	srand(time(NULL));
	/* create entire population */
	for (i=0; i<size; ++i) {
		/* create each genome */
		pop->pool[i].wts = (char *)calloc(gene_cnt, sizeof(char));
		for (j=0; j<gene_cnt; ++j)
			/* initialize each genome with random gene values */
			pop->pool[i].wts[j] = rand() % (UPPER - LOWER + 1) + LOWER;
	}
	return pop;
}

/* destroy genetic algorithm and deallocate memory */
void ga_destroy(pop_t *const pop)
{
	int i;
	for (i=0; i< pop->size; ++i)
		free(pop->pool[i].wts);					/* free each genome */
	free(pop->pool);							/* free population pool */
	free(pop->pwt1);							/* free parent 1 genes */
	free(pop->pwt2);							/* free parent 2 genes */
	free(pop->cwt1);							/* free child 1 genes */
	free(pop->cwt2);							/* free child 2 genes */
	free(pop);
}

/* one population life cycle */
void ga_epoch(pop_t *const pop)
{
	char *wts[pop->size];
	int i, j;
	for (i=0; i< pop->size; ++i)
		wts[i] = (char *)calloc(pop->size, sizeof(char));
	pop_error(pop);
	qsort(pop->pool, pop->size, sizeof(genome_t), cmp_genome);
	for (i=0; i< pop->size; i+=2) {
		tourn_select(pop, pop->pwt1, pop->pwt2);
		crossover(pop, pop->pwt1, pop->pwt2, pop->cwt1, pop->cwt2);
		mutate(pop, pop->cwt1);
		mutate(pop, pop->cwt2);
		for (j=0; j< pop->gene_cnt; ++j) {
			wts[i][j] = pop->cwt1[j];
			wts[i+1][j] = pop->cwt2[j];
		}
	}
	for (i=0; i< pop->size; ++i) {
		for (j=0; j< pop->gene_cnt; ++j)
			pop->pool[i].wts[j] = wts[i][j];
		free(wts[i]);
	}
}

/* run genetic algorithm for several generations */
genome_t *ga_run(pop_t *const pop, const int num_gen)
{
	int i;
	for (i=0; i<num_gen; ++i) {
		ga_epoch(pop);
		printf("\rbest:   [ %s ] error: %d generation: %d",
			pop->best->wts, pop->best->err, i+1);
		if (pop->found == true) {				/* termination condition */
			return pop->best;
		}
	}
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

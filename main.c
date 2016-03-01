#include "main.h"

int main(void)
{
	printf("[ genetic algorithm ]\n");

	/* target string */
	const char target[] = "Shraman Biswas";

	unsigned int num_gen, gene_cnt, size, tourn_size;
	double c_rate, m_rate;
	genome_t *best=NULL;
	pop_t *pop=NULL;

	/* genetic algorithm parameters */
	num_gen = 100000;			/* number of generations */
	gene_cnt = SIZE(target) - 1;		/* number of genes */
	size = 20;				/* population size */
	c_rate = 0.9;				/* crossover rate */
	m_rate = 0.1;				/* mutation rate */
	tourn_size = 3;				/* tournament size */

	/* display target string */
	printf("target: [ %s ]\n", target);

	/* create genetic algorithm and obtain population */
	pop = ga_create(target, gene_cnt, size, c_rate, m_rate, tourn_size);

	/* apply genetic algorithm and obtain best genome if found */
	best = ga_run(pop, num_gen);

	/* check if best genome was found */
	if (best != NULL) {
		printf("\ngenome found!\n");
	} else {
		printf("\ngenome not found!\n")	;
	}

	/* destroy genetic algorithm population */
	ga_destroy(pop);

	return EXIT_SUCCESS;
}

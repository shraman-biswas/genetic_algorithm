#include "main.h"

int main()
{
	printf("[ genetic algorithm ]\n");

	const char target[] = "SHRAMAN";

	bool found=false;
	int num_gen, gene_cnt, size, tourn_size;
	double c_rate, m_rate;
	pop_t *pop=NULL;

	/* genetic algorithm parameters */
	num_gen = 10000;							/* number of generations */
	gene_cnt = SIZE(target) - 1;				/* number of genes */
	size = 10;									/* population size */
	c_rate = 0.9;								/* crossover rate */
	m_rate = 0.1;								/* mutation rate */
	tourn_size = 3;								/* tournament size */

	/* create genetic algorithm and obtain population */
	pop = ga_create(target, gene_cnt, size, c_rate, m_rate, tourn_size);

	/* apply genetic algorithm */
	found = ga_run(pop, num_gen);

	/* check if genome was found */
	if (found == false)
		printf("\ngenome not found!\n");

	/* destroy genetic algorithm population */
	ga_destroy(pop);

	return EXIT_SUCCESS;
}

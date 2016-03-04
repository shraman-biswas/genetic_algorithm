Genetic Algorithm
=================

Genetic algorithm implementation with tournament selection.  
Performs basic string matching. Evolves the fittest genome string that matches the target string perfectly.  
Uses simple crossover operation that swaps parent genes at a randomly generated crossover point.  
Uses simple mutation operation that increments or decrements every character in genome string.   
Uses tournament selection to select 2 fit parents from the entire population pool.

Features:
* genetic algorithm creation
* genetic algorithm destruction
* genetic algorithm epoch
* genetic algorithm run
* display genome helper function

Parameters:
* target string
* number of genes in each genome (length of target string)
* number of genomes (population size)
* crossover rate
* mutation rate
* tournament selection size

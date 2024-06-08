/*
 * Exact genetic sequence alignment
 * (Using brute force)
 *
 * OpenMP version
 *
 * Computacion Paralela, Grado en Informatica (Universidad de Valladolid)
 * 2023/2024
 *
 * v1.0.1
 *
 * (c) 2024, Arturo Gonzalez-Escribano
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<sys/time.h>
#include<omp.h>


/* Arbitrary value to indicate that no matches are found */
#define	NOT_FOUND	-1

/* Arbitrary value to restrict the checksums period */
#define CHECKSUM_MAX	65536

/*
 * Function: Increment the number of pattern matches on the sequence positions
 * 	This function can be changed and/or optimized by the students
 */
/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */
void increment_matches( int pat, unsigned long *pat_found, long *pat_length, int *seq_matches ) {
	
	/*AL INTENTAR PARALELIZAR ESTE BUCLE EL TIEMPO ME DA PEOR QUE SIN PARALELIZAR :/ */
	//He probado con ifs dependiendo del pat_length y tampoco mejora
	//creo que con lo que pone el codigo de arriba se tiene que cambiar pero no optimizar con los pragmas...
	//Si sacamos las variables fuera las comprueba una sola vez, no cada vez que recorre el bucle
	int offset = pat_found[pat];
	int limit = pat_length[pat];
	int ind;

	for (int ind = 0; ind < limit; ind++) {
		if (seq_matches[offset + ind] == NOT_FOUND) {
			seq_matches[offset + ind] = 0;
		} else {
			seq_matches[offset + ind]++;
		}
	}
}
/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */


/* 
 * Utils: Function to get wall time
 */
double cp_Wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.0e-6 * tv.tv_usec;
}


/*
 * Utils: Random generator
 */
#include "rng.c"


/*
 * Function: Allocate new patttern
 */
char *pattern_allocate( rng_t *random, unsigned long pat_rng_length_mean, unsigned long pat_rng_length_dev, unsigned long seq_length, unsigned long *new_length ) {

	/* Random length */
	double length = (unsigned long)rng_next_normal( random, (double)pat_rng_length_mean, (double)pat_rng_length_dev );
	if ( length > seq_length ) length = seq_length;
	if ( length <= 0 ) length = 1;

	/* Allocate pattern */
	char *pattern = (char *)malloc( sizeof(char) * length );
	if ( pattern == NULL ) {
		fprintf(stderr,"\n-- Error allocating a pattern of size: %lu\n", length );
		exit( EXIT_FAILURE );
	}

	/* Return results */
	*new_length = length;
	return pattern;
}

/*
 * Function: Fill random sequence or pattern
 */
void generate_rng_sequence( rng_t *random, float prob_G, float prob_C, float prob_A, char *seq, unsigned long length) {
	#pragma omp parallel
	{
		rng_t local_random = *random;
		int first_iteration = 1;

		unsigned long ind; 
		#pragma omp for
		for( ind=0; ind<length; ind++ ) {
			if ( first_iteration ) {
				rng_skip( &local_random, ind );
				first_iteration = 0;
			}
			double prob = rng_next( &local_random );
			if( prob < prob_G ) seq[ind] = 'G';
			else if( prob < prob_C ) seq[ind] = 'C';
			else if( prob < prob_A ) seq[ind] = 'A';
			else seq[ind] = 'T';
		}
	}
	rng_skip( random, length );
}

/*
 * Function: Copy a sample of the sequence
 */
void generate_sample_sequence( rng_t *random, char *sequence, unsigned long seq_length, unsigned long pat_samp_loc_mean, unsigned long pat_samp_loc_dev, char *pattern, unsigned long length) {
	/* Choose location */
	unsigned long  location = (unsigned long)rng_next_normal( random, (double)pat_samp_loc_mean, (double)pat_samp_loc_dev );
	if ( location > seq_length - length ) location = seq_length - length;
	if ( location <= 0 ) location = 0;

	/* Copy sample */
	unsigned long ind; 
	#pragma omp parallel for
	for( ind=0; ind<length; ind++ )
		pattern[ind] = sequence[ind+location];
}


/*
 * Function: Print usage line in stderr
 */
void show_usage( char *program_name ) {
	fprintf(stderr,"Usage: %s ", program_name );
	fprintf(stderr,"<seq_length> <prob_G> <prob_C> <prob_A> <pat_rng_num> <pat_rng_length_mean> <pat_rng_length_dev> <pat_samples_num> <pat_samp_length_mean> <pat_samp_length_dev> <pat_samp_loc_mean> <pat_samp_loc_dev> <pat_samp_mix:B[efore]|A[fter]|M[ixed]> <long_seed>\n");
	fprintf(stderr,"\n");
}



/*
 * MAIN PROGRAM
 */
int main(int argc, char *argv[]) {
	/* 0. Default output and error without buffering, forces to write immediately */
	setbuf(stdout, NULL);
	setbuf(stderr, NULL);

	/* 1. Read scenary arguments */
	/* 1.1. Check minimum number of arguments */
	if (argc < 15) {
		fprintf(stderr, "\n-- Error: Not enough arguments when reading configuration from the command line\n\n");
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}

	/* 1.2. Read argument values */
	unsigned long seq_length = atol( argv[1] );
	float prob_G = atof( argv[2] );
	float prob_C = atof( argv[3] );
	float prob_A = atof( argv[4] );
	if ( prob_G + prob_C + prob_A > 1 ) {
		fprintf(stderr, "\n-- Error: The sum of G,C,A,T nucleotid probabilities cannot be higher than 1\n\n");
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}
	prob_C += prob_G;
	prob_A += prob_C;

	int pat_rng_num = atoi( argv[5] );
	unsigned long pat_rng_length_mean = atol( argv[6] );
	unsigned long pat_rng_length_dev = atol( argv[7] );
	
	int pat_samp_num = atoi( argv[8] );
	unsigned long pat_samp_length_mean = atol( argv[9] );
	unsigned long pat_samp_length_dev = atol( argv[10] );
	unsigned long pat_samp_loc_mean = atol( argv[11] );
	unsigned long pat_samp_loc_dev = atol( argv[12] );

	char pat_samp_mix = argv[13][0];
	if ( pat_samp_mix != 'B' && pat_samp_mix != 'A' && pat_samp_mix != 'M' ) {
		fprintf(stderr, "\n-- Error: Incorrect first character of pat_samp_mix: %c\n\n", pat_samp_mix);
		show_usage( argv[0] );
		exit( EXIT_FAILURE );
	}

	unsigned long seed = atol( argv[14] );

#ifdef DEBUG
	/* DEBUG: Print arguments */
	printf("\nArguments: seq_length=%lu\n", seq_length );
	printf("Arguments: Accumulated probabilitiy G=%f, C=%f, A=%f, T=1\n", prob_G, prob_C, prob_A );
	printf("Arguments: Random patterns number=%d, length_mean=%lu, length_dev=%lu\n", pat_rng_num, pat_rng_length_mean, pat_rng_length_dev );
	printf("Arguments: Sample patterns number=%d, length_mean=%lu, length_dev=%lu, loc_mean=%lu, loc_dev=%lu\n", pat_samp_num, pat_samp_length_mean, pat_samp_length_dev, pat_samp_loc_mean, pat_samp_loc_dev );
	printf("Arguments: Type of mix: %c, Random seed: %lu\n", pat_samp_mix, seed );
	printf("\n");
#endif // DEBUG

	/* 2. Initialize data structures */
	/* 2.1. Allocate and fill sequence */
	char *sequence = (char *)malloc( sizeof(char) * seq_length );
	if ( sequence == NULL ) {
		fprintf(stderr,"\n-- Error allocating the sequence for size: %lu\n", seq_length );
		exit( EXIT_FAILURE );
	}
	rng_t random = rng_new( seed );
	generate_rng_sequence( &random, prob_G, prob_C, prob_A, sequence, seq_length);

	/* 2.2. Allocate and fill patterns */
	/* 2.2.1 Allocate main structures */
	int pat_number = pat_rng_num + pat_samp_num;
	char **patterns;
	unsigned long *pat_length = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	char **pattern = (char **)malloc( sizeof(char*) * pat_number );
	if ( pattern == NULL || pat_length == NULL ) {
		fprintf(stderr,"\n-- Error allocating the basic patterns structures for size: %d\n", pat_number );
		exit( EXIT_FAILURE );
	}

	/* 2.2.2 Allocate and initialize ancillary structure for pattern types */
	unsigned long ind;
	#define PAT_TYPE_NONE	0
	#define PAT_TYPE_RNG	1
	#define PAT_TYPE_SAMP	2
	char *pat_type = (char *)malloc( sizeof(char) * pat_number );
	if ( pat_type == NULL ) {
		fprintf(stderr,"\n-- Error allocating ancillary structure for pattern of size: %d\n", pat_number );
		exit( EXIT_FAILURE );
	}
	for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_NONE;

	/* 2.2.3 Fill up pattern types using the chosen mode */
	switch( pat_samp_mix ) {
	case 'A':
		for( ind=0; ind<pat_rng_num; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		for( ; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		break;
	case 'B':
		for( ind=0; ind<pat_samp_num; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		for( ; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		break;
	default:
		if ( pat_rng_num == 0 ) {
			for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_SAMP;
		}
		else if ( pat_samp_num == 0 ) {
			for( ind=0; ind<pat_number; ind++ ) pat_type[ind] = PAT_TYPE_RNG;
		}
		else if ( pat_rng_num < pat_samp_num ) {
			int interval = pat_number / pat_rng_num;
			for( ind=0; ind<pat_number; ind++ ) 
				if ( (ind+1) % interval == 0 ) pat_type[ind] = PAT_TYPE_RNG;
				else pat_type[ind] = PAT_TYPE_SAMP;
		}
		else {
			int interval = pat_number / pat_samp_num;
			for( ind=0; ind<pat_number; ind++ ) 
				if ( (ind+1) % interval == 0 ) pat_type[ind] = PAT_TYPE_SAMP;
				else pat_type[ind] = PAT_TYPE_RNG;
		}
	}

	/* 2.2.4 Generate the patterns */
	for( ind=0; ind<pat_number; ind++ ) {
		if ( pat_type[ind] == PAT_TYPE_RNG ) {
			pattern[ind] = pattern_allocate( &random, pat_rng_length_mean, pat_rng_length_dev, seq_length, &pat_length[ind] );
			generate_rng_sequence( &random, prob_G, prob_C, prob_A, pattern[ind], pat_length[ind] );
		}
		else if ( pat_type[ind] == PAT_TYPE_SAMP ) {
			pattern[ind] = pattern_allocate( &random, pat_samp_length_mean, pat_samp_length_dev, seq_length, &pat_length[ind] );
			generate_sample_sequence( &random, sequence, seq_length, pat_samp_loc_mean, pat_samp_loc_dev, pattern[ind], pat_length[ind] );
		}
		else {
			fprintf(stderr,"\n-- Error internal: Paranoic check! A pattern without type at position %d\n", ind );
			exit( EXIT_FAILURE );
		}
	}
	free( pat_type );

#ifdef DEBUG
	/* DEBUG: Print sequence and patterns */
	printf("-----------------\n");
	printf("Sequence: ");
	for( ind=0; ind<seq_length; ind++ ) 
		printf( "%c", sequence[ind] );
	printf("\n-----------------\n");
	printf("Patterns: %d ( rng: %d, samples: %d )\n", pat_number, pat_rng_num, pat_samp_num );
	int debug_pat;
	for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
		printf( "Pat[%lu]: ", debug_pat );
		for( ind=0; ind<pat_length[debug_pat]; ind++ ) 
			printf( "%c", pattern[debug_pat][ind] );
		printf("\n");
	}
	printf("-----------------\n\n");
#endif // DEBUG

	/* Avoid the usage of arguments to take strategic decisions
	 * In a real case the user only has the patterns and sequence data to analize
	 */
	argc = 0;
	argv = NULL;
	prob_G = 0.0;
	prob_C = 0.0;
	prob_A = 0.0;
	pat_rng_num = 0;
	pat_rng_length_mean = 0;
	pat_rng_length_dev = 0;
	pat_samp_num = 0;
	pat_samp_length_mean = 0;
	pat_samp_length_dev = 0;
	pat_samp_loc_mean = 0;
	pat_samp_loc_dev = 0;
	seed = 0;

	/* 2.3. Other result data and structures */
	int pat_matches = 0;
	unsigned long *pat_found;
	int *seq_matches;
	int *seq_longest;

	pat_found = (unsigned long *)malloc( sizeof(unsigned long) * pat_number );
	seq_matches = (int *)malloc( sizeof(int) * seq_length );
	seq_longest = (int *)malloc( sizeof(int) * seq_length );
	if ( seq_matches == NULL || seq_longest == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux sequence structures for size: %lu\n", seq_length );
		exit( EXIT_FAILURE );
	}
	if ( pat_found == NULL ) {
		fprintf(stderr,"\n-- Error allocating aux pattern structure for size: %lu\n", pat_length );
		exit( EXIT_FAILURE );
	}

	
	/* 3. Start global timer */
	double ttotal = cp_Wtime();

/*
 *
 * START HERE: DO NOT CHANGE THE CODE ABOVE THIS POINT
 *
 */
	

	/* 4. Initialize ancillary structures */
	/* ESTOS DOS FOR ESTAN BASTANTE CERCA Y HACEN CASI LO MISMO POR LO QUE 
	   HE HECHO UN REGION PARALEA Y 2 FOR PARA UTIL LOS MISMOS THREADS*/
	
	#pragma omp parallel shared(pat_found,seq_matches,seq_longest) private(ind) 
	{	
		#pragma omp for schedule(static)
		for( ind=0; ind<pat_number; ind++) {
			pat_found[ind] = NOT_FOUND;
		}
		#pragma omp for schedule(static)
		for( ind=0; ind<seq_length; ind++) {
			seq_matches[ind] = 0;
			seq_longest[ind] = 0;
		}
	}
	

	/* 5. Search for each pattern */


	unsigned long start;
	unsigned long pat;
	#pragma omp parallel for shared(pat_length,sequence,pattern,pat_found) private(pat,start,ind) schedule(static)
	for( pat=0; pat < pat_number; pat++ ) {
		/* 5.1. For each posible starting position */
		unsigned long pl = pat_length[pat];
		unsigned long slpl = seq_length - pl;

		for( start=0; start <= slpl; start++) {
			/* 5.1.1. For each pattern element */
			for( ind=0; ind<pl; ind++) {
				/* Stop this test when different nucleotids are found */
				if ( sequence[start + ind] != pattern[pat][ind] ) break;
			}
			/* 5.1.2. Check if the loop ended with a match */
			if ( ind == pl ) {
				pat_found[pat] = start;
				break;
			
		}

		/* 5.2. Pattern found */
		if ( pat_found[pat] != -1 ) {
			/* 4.2.1. Increment the number of pattern matches on the sequence positions */
			#pragma omp critical
			{
				pat_matches++;
				int limit = pat_length[pat];
				int ind;

				for (int ind = 0; ind < limit; ind++) {
					if (seq_matches[start + ind] != -1) {
						seq_matches[start + ind]++;
					} else {
						seq_matches[start + ind] = 0;
					}
				}
			}
		}
	}

	

	/* 6. Annotate the index of the longest pattern matched on each position */
	/* UN FOR PARALELO NORMAL, PARECE TODO PARALIZABLE*/
	//Saco el rellenado del vector fuera para alomejor poder colapsar

		
	for( pat=0; pat<pat_number; pat++ )  {
		
		#pragma omp task
		{
			int p_f = pat_found[pat];
			if ( p_f != -1 ){				
				#pragma omp parallel for shared(seq_longest,pat_found,pat_length) private(ind) schedule(static)
					for( ind=0; ind < seq_length; ind++){
							if ( p_f <= ind){
								int p_l = pat_length[pat];
								int p_t = p_f + p_l;
								if(ind < p_t ){
									int s_l = seq_longest[ind] ;
									if ( s_l < p_l )
										seq_longest[ind] = pat_length[pat];
								}
							}
					}
				}	
		}
	}

	

	/* 7. Check sums */
	/*SE HA HECHO LA REDUCCION PARA TODAS LAS VARIABLES*/

		

	unsigned long checksum_matches = 0;
	unsigned long checksum_longest = 0;
	unsigned long checksum_found = 0;

	
	#pragma omp parallel for shared(pat_found) reduction(+:checksum_found) schedule(static)
	for( ind=0; ind < pat_number; ind++) {
		if ( pat_found[ind] != -1 )
			checksum_found = ( checksum_found + pat_found[ind] );
	}
	checksum_found = checksum_found % 65536;

	#pragma omp parallel for shared(seq_matches) reduction(+:checksum_matches) schedule(static)
	for( ind=0; ind < seq_length; ind++) {
		if ( seq_matches[ind] != -1 )
			checksum_matches = ( checksum_matches + seq_matches[ind] );
	}
	checksum_matches = checksum_matches % 65536;
	
	#pragma omp parallel for shared(seq_longest) reduction(+:checksum_longest) schedule(static)
	for( ind=0; ind < seq_length; ind++) {
		checksum_longest = ( checksum_longest + seq_longest[ind] );
	}
	checksum_longest = checksum_longest % 65536;

/*
 *
 * STOP HERE: DO NOT CHANGE THE CODE BELOW THIS POINT
 *
 */

#ifdef DEBUG
	/* DEBUG: Write results */
	printf("-----------------\n");
	printf("Found start:");
	for( debug_pat=0; debug_pat<pat_number; debug_pat++ ) {
		printf( " %lu", pat_found[debug_pat] );
	}
	printf("\n");
	printf("-----------------\n");
	printf("Matches:");
	for( ind=0; ind<seq_length; ind++ ) 
		printf( " %d", seq_matches[ind] );
	printf("\n");
	printf("-----------------\n");
	printf("Longest:");
	for( ind=0; ind<seq_length; ind++ ) 
		printf( " %d", seq_longest[ind] );
	printf("\n");
	printf("-----------------\n");
#endif // DEBUG

	/* 8. Stop global time */
	ttotal = cp_Wtime() - ttotal;

	/* 9. Output for leaderboard */
	printf("\n");
	/* 9.1. Total computation time */
	printf("Time: %lf\n", ttotal);

	/* 9.2. Results: Statistics */
	printf("Result: %d, %lu, %lu, %lu\n\n", 
			pat_matches,
			checksum_found,
			checksum_matches,
			checksum_longest );
		
	/* 10. Free resources */	
	int i;
	for( i=0; i<pat_number; i++ ) free( pattern[i] );
	free( pattern );
	free( sequence );
	free( pat_length );
	free( pat_found );
	free( seq_longest );
	free( seq_matches );

	/* 11. End */
	return 0;
}

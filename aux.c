	double startTime = omp_get_wtime();

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
	
	double endTime = omp_get_wtime();

	/* 5. Search for each pattern */
	/* 5. Search for each pattern */

	unsigned long start;
	unsigned long pat;
	startTime = omp_get_wtime();
	#pragma omp parallel for shared(pat_length,sequence,pattern,pat_found) private(pat,start,ind) schedule(static)
	for( pat=0; pat < pat_number; pat++ ) {
		/* 5.1. For each posible starting position */
		for( start=0; start <= seq_length - pat_length[pat]; start++) {

			/* 5.1.1. For each pattern element */
			for( ind=0; ind<pat_length[pat]; ind++) {
				/* Stop this test when different nucleotids are found */
				if ( sequence[start + ind] != pattern[pat][ind] ) break;
			}
			/* 5.1.2. Check if the loop ended with a match */
			if ( ind == pat_length[pat] ) {
				#pragma omp atomic
				pat_matches++;
				pat_found[pat] = start;
				break;
			}
		}

		/* 5.2. Pattern found */
		if ( pat_found[pat] != NOT_FOUND ) {
			/* 4.2.1. Increment the number of pattern matches on the sequence positions */
			#pragma omp critical
			increment_matches( pat, pat_found, pat_length, seq_matches );
		}
	}

	
	endTime = omp_get_wtime();

	/* 6. Annotate the index of the longest pattern matched on each position */
	/* UN FOR PARALELO NORMAL, PARECE TODO PARALIZABLE*/
	//Saco el rellenado del vector fuera para alomejor poder colapsar

		
	startTime = omp_get_wtime();
	memset(seq_longest, 0, sizeof(seq_longest));
	#pragma omp parallel for shared(seq_longest,pat_found,pat_length) private(ind,pat) schedule(dynamic,12)
	for( ind=0; ind < seq_length; ind++) {

		for( pat=0; pat<pat_number; pat++ ) {
			
			int p_f = pat_found[pat];
			int p_l = pat_length[pat];
			int p_t = p_f + p_l;
			int s_l = seq_longest[ind] ;
			if ( p_f != NOT_FOUND )
				if ( p_f <= ind)
					if(ind < p_t )
				//Alomejor aqui se puede usar una reduction(min:seq_longest[ind])
						if ( s_l < p_l )
							seq_longest[ind] = pat_length[pat];
			}
		
	}

		
	endTime = omp_get_wtime();

	/* 7. Check sums */
	/*SE HA HECHO LA REDUCCION PARA TODAS LAS VARIABLES*/

		
	startTime = omp_get_wtime();

	unsigned long checksum_matches = 0;
	unsigned long checksum_longest = 0;
	unsigned long checksum_found = 0;

	#pragma omp parallel for shared(pat_found) reduction(+:checksum_found) schedule(static)
	for( ind=0; ind < pat_number; ind++) {
		if ( pat_found[ind] != NOT_FOUND )
			checksum_found = ( checksum_found + pat_found[ind] );
	}
	checksum_found = checksum_found % CHECKSUM_MAX;

	#pragma omp parallel for shared(seq_matches) reduction(+:checksum_matches) schedule(static)
	for( ind=0; ind < seq_length; ind++) {
		if ( seq_matches[ind] != NOT_FOUND )
			checksum_matches = ( checksum_matches + seq_matches[ind] );
	}
	checksum_matches = checksum_matches % CHECKSUM_MAX;
	
	#pragma omp parallel for shared(seq_longest) reduction(+:checksum_longest) schedule(static)
	for( ind=0; ind < seq_length; ind++) {
		checksum_longest = ( checksum_longest + seq_longest[ind] );
	}
	checksum_longest = checksum_longest % CHECKSUM_MAX;

	endTime = omp_get_wtime();

#include "edwards_curve.h"

void random_key(uint8_t key[])
{
	uint8_t i, tmp, r;
	int8_t exp, sgn;
	for(i = 0; i < N; i++)
	{

		r = B[i] & 0x1;						// B_i mod 2

		// exp is randomly selected from |[ 0, B ]|
		randombytes(&tmp, 1);
		while ( issmaller((int32_t)B[i], (int32_t)tmp) == -1 )	// constant-time comparison
			randombytes(&tmp, 1);

		exp = (int8_t)tmp;

		// Mapping integers from |[ 0, B |] into
		//                                      |[ -B/2, B/2]| if B is even, or
		//                                      |[ -(B+1)/2, (B-1)/2 ]| if B is odd.
		exp = ( (exp << 1) - (B[i] + r) ) >> 1;

		// Mapping into the set |[-B, B]|.
		exp = (exp << 1) + r;
		sgn = exp >> 7;	// sign of exp

		// Next, to write  key[i] = e || ((1 + sgn)/2)
		cmov(&exp, -exp, sgn == -1);
		key[i] = (exp << 1) ^ (1 & (1 + sgn));
	};
};

void printf_key(uint8_t key[], char *c)
{
	int i;
	printf("%s := ", c);
	printf("{\t  %3d", (int)( (2*(key[0] & 0x1) - 1) * (key[0] >> 1) ));

	for(i = 1; i < N; i++)
	{
		printf(", %3d", (int)( (2*(key[i] & 0x1) - 1) * (key[i] >> 1) ) );
		if( (i % 18) == 17 )
			printf("\n\t\t");
	};

	printf("};\n");

};

/* ------------------------------------------------------------------------------------------------------------------- *
   strategy_evaluation()
   inputs: a secret key, the projective Edwards y-coordinates of y(T_{+}) = YT_{+}/ZT_{+} and y(T_{-}) = YT_{-}/ZT_{-}, 
           the Edwards curve constants A[0]:=a, and A[1]:=(a - d), and an integer number 0 <= i <= ROUNDS;
   output: the isogenous Edwards curve constants B[0]:=a and B[1]:=(a - d) determined by the i-th strategy evaluated at 
           segment of the secret key
   
   NOTE: the point T_{+} must has order less than or equal to (p+1)/4, and this is the algorithm that uses dummy operations
 * ------------------------------------------------------------------------------------------------------------------- */
void strategy_evaluation(proj C, uint8_t counter[], uint8_t key[], const proj T_p, const proj T_m, const proj A, const int j)
{
	proj current_A, 				// Current value of a and (a -d)
	     tmp_Ts[(L[L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1]] >> 1)];	// This will be used to evaluate points under isogenies
	point_copy(current_A, A);			// a and (a - d)

	proj SPLITTING_POINTS[2][LOG2_OF_N_PLUS_ONE];	// This will correspond with the splitting points of the strategy
	proj G[2];

	point_copy(SPLITTING_POINTS[0][0], T_m);	// T_{-}
	point_copy(SPLITTING_POINTS[1][0], T_p);	// T_{+}

	int STRATEGY_ELEMENT = 0, 			// Current element of the strategy to be used
	    local_i, local_j;

	int BLOCK = 0,					// BLOCK is used for determined when a point has order l
	    SPLIT_ELEMENT = 0;				// At the beginning we have only one point in each split

	int yMUL_PERFORMED[LOG2_OF_N_PLUS_ONE];		// The current number of scalar multiplication performed

	uint32_t bc;
	int8_t ec = 0, mask = 0;

	// ------------------------------------------------------------------------------------------------------------------------
	// Main loop: constructing degree-l_j isogenies with k = 1, 2, ..., NUMBER_OF_PRIMES[j] - 1
	for(local_j = 0; local_j < (NUMBER_OF_PRIMES[j] - 1); local_j++)
	{
        ec = lookup(L_STRATEGY[j][local_j], key);	// To get current e_i in constant-time
		// Next, the split of points are computed and saved them in order to be evaluated and used in the next isogeny construction
        while( (BLOCK + S[j][STRATEGY_ELEMENT]) < (NUMBER_OF_PRIMES[j] -  1 - local_j) )
        {

            // A split will be added
            SPLIT_ELEMENT += 1;

            // We set the seed of the new split to be computed and saved
            point_copy(SPLITTING_POINTS[0][SPLIT_ELEMENT], SPLITTING_POINTS[0][SPLIT_ELEMENT - 1]);
            point_copy(SPLITTING_POINTS[1][SPLIT_ELEMENT], SPLITTING_POINTS[1][SPLIT_ELEMENT - 1]);
            for(local_i = BLOCK; local_i < (BLOCK + S[j][STRATEGY_ELEMENT]); local_i++)
            {
                // Split corresponding to T_{-}
                yMUL(SPLITTING_POINTS[0][SPLIT_ELEMENT], SPLITTING_POINTS[0][SPLIT_ELEMENT], current_A, L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1 - local_i]);
                // Split corresponding to T_{+}
                yMUL(SPLITTING_POINTS[1][SPLIT_ELEMENT], SPLITTING_POINTS[1][SPLIT_ELEMENT], current_A, L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1 - local_i]);
            };
            
            yMUL_PERFORMED[SPLIT_ELEMENT] = S[j][STRATEGY_ELEMENT];	// The number of yMUL performed is saved
            BLOCK += S[j][STRATEGY_ELEMENT];					// BLOCK is increased by the number of yMUL performed
            STRATEGY_ELEMENT += 1;							// Next, we move to the next element of the strategy
        }

        if( BLOCK < (NUMBER_OF_PRIMES[j] -  1 - local_j) )
        {
            SPLIT_ELEMENT += 1;
            point_copy(SPLITTING_POINTS[0][SPLIT_ELEMENT], SPLITTING_POINTS[0][SPLIT_ELEMENT - 1]);
            point_copy(SPLITTING_POINTS[1][SPLIT_ELEMENT], SPLITTING_POINTS[1][SPLIT_ELEMENT - 1]);
                
            fp_cswap(SPLITTING_POINTS[0][SPLIT_ELEMENT][0], SPLITTING_POINTS[1][SPLIT_ELEMENT][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
            fp_cswap(SPLITTING_POINTS[0][SPLIT_ELEMENT][1], SPLITTING_POINTS[1][SPLIT_ELEMENT][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
            
            for(local_i = BLOCK; local_i < (BLOCK + S[j][STRATEGY_ELEMENT]); local_i++)
                yMUL(SPLITTING_POINTS[0][SPLIT_ELEMENT], SPLITTING_POINTS[0][SPLIT_ELEMENT], current_A, L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1 - local_i]);
            
            fp_cswap(SPLITTING_POINTS[0][SPLIT_ELEMENT][0], SPLITTING_POINTS[1][SPLIT_ELEMENT][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
            fp_cswap(SPLITTING_POINTS[0][SPLIT_ELEMENT][1], SPLITTING_POINTS[1][SPLIT_ELEMENT][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
            
            yMUL_PERFORMED[SPLIT_ELEMENT] = S[j][STRATEGY_ELEMENT];	// The number of yMUL performed is saved
            BLOCK += S[j][STRATEGY_ELEMENT];					// BLOCK is increased by the number of yMUL performed
            STRATEGY_ELEMENT += 1;							// Next, we move to the next element of the strategy
        }
		
        // Now, a degree-(l_{batches[m][i]}) will be constructed. Let l = l_{batches[m][i]}.
		point_copy(G[0], SPLITTING_POINTS[0][SPLIT_ELEMENT]);	// order-l point determined by T_{-}
		point_copy(G[1], SPLITTING_POINTS[1][SPLIT_ELEMENT]);	// order-l point determined by T_{+}
        
        fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
        fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

        if(isequal(counter[L_STRATEGY[j][local_j]], B[L_STRATEGY[j][local_j]]) != 1)
        {
    		if(isinfinity(G[0]) != 1)
    		{
    			// At this point, BLOCK = NUMBER_OF_PRIMES[j] - 1 - local_j implies that the last two points in the split have order at most l_{local_j}
    			bc = isequal(ec >> 1, 0) & 1;			// Bit that determines if a dummy operation will be perfomed
    			key[L_STRATEGY[j][local_j]] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
    			counter[L_STRATEGY[j][local_j]] += 1;

    			// At this point, our kernel has order l_{local_j}. Therefore, we can construct a degree-l_{local_j} isogeny
    			yISOG(tmp_Ts, current_A, SPLITTING_POINTS[(ec & 0x1)][SPLIT_ELEMENT], current_A, L_STRATEGY[j][local_j]);
    			for(local_i = 0; local_i < SPLIT_ELEMENT; local_i++)
    			{

                    fp_cswap(SPLITTING_POINTS[0][local_i][0], SPLITTING_POINTS[1][local_i][0], (ec & 1));
                    fp_cswap(SPLITTING_POINTS[0][local_i][1], SPLITTING_POINTS[1][local_i][1], (ec & 1));
                    
    				// Evaluating the split corresponding with T_{-}
    				yEVAL(SPLITTING_POINTS[0][local_i], SPLITTING_POINTS[0][local_i], tmp_Ts, L_STRATEGY[j][local_j]);
    				// Evaluating the split corresponding with T_{+}
    				yEVAL(SPLITTING_POINTS[1][local_i], SPLITTING_POINTS[1][local_i], tmp_Ts, L_STRATEGY[j][local_j]);
    				// In additon, the split that doens't correspond with the kernel must be multiplied by [l_{local_j}]
    				yMUL(SPLITTING_POINTS[1][local_i], SPLITTING_POINTS[1][local_i], current_A, L_STRATEGY[j][local_j]);
                                    
                    fp_cswap(SPLITTING_POINTS[0][local_i][0], SPLITTING_POINTS[1][local_i][0], (ec & 1));
                    fp_cswap(SPLITTING_POINTS[0][local_i][1], SPLITTING_POINTS[1][local_i][1], (ec & 1));

    			};
    		}
    		else
    		{
    			for(local_i = 0; local_i < SPLIT_ELEMENT; local_i++)
    			{
                    fp_cswap(SPLITTING_POINTS[0][local_i][0], SPLITTING_POINTS[1][local_i][0], (ec & 1));
                    fp_cswap(SPLITTING_POINTS[0][local_i][1], SPLITTING_POINTS[1][local_i][1], (ec & 1));
                    
    				// In additon, the split that doens't correspond with the kernel must be multiplied by [l_{local_j}]
    				yMUL(SPLITTING_POINTS[1][local_i], SPLITTING_POINTS[1][local_i], current_A, L_STRATEGY[j][local_j]);
                    
                    fp_cswap(SPLITTING_POINTS[0][local_i][0], SPLITTING_POINTS[1][local_i][0], (ec & 1));
                    fp_cswap(SPLITTING_POINTS[0][local_i][1], SPLITTING_POINTS[1][local_i][1], (ec & 1));
    			};
    		};
        }
        else
        {
            // This branch only depends on randomness
            if(isinfinity(G[0]) != 1)
            {
                for(local_i = 0; local_i < SPLIT_ELEMENT; local_i++)
                {
                    yMUL(SPLITTING_POINTS[0][local_i], SPLITTING_POINTS[0][local_i], current_A, L_STRATEGY[j][local_j]);
                    yMUL(SPLITTING_POINTS[1][local_i], SPLITTING_POINTS[1][local_i], current_A, L_STRATEGY[j][local_j]);
                };
            }
            else
            {
    			for(local_i = 0; local_i < SPLIT_ELEMENT; local_i++)
    			{
                    fp_cswap(SPLITTING_POINTS[0][local_i][0], SPLITTING_POINTS[1][local_i][0], (ec & 1));
                    fp_cswap(SPLITTING_POINTS[0][local_i][1], SPLITTING_POINTS[1][local_i][1], (ec & 1));
                    
    				// In additon, the split that doens't correspond with the kernel must be multiplied by [l_{local_j}]
    				yMUL(SPLITTING_POINTS[1][local_i], SPLITTING_POINTS[1][local_i], current_A, L_STRATEGY[j][local_j]);
                    
                    fp_cswap(SPLITTING_POINTS[0][local_i][0], SPLITTING_POINTS[1][local_i][0], (ec & 1));
                    fp_cswap(SPLITTING_POINTS[0][local_i][1], SPLITTING_POINTS[1][local_i][1], (ec & 1));
    			};
            }
        };

		BLOCK -= yMUL_PERFORMED[SPLIT_ELEMENT];				// BLOCK is decreased by the last number of yMUL performed
		yMUL_PERFORMED[SPLIT_ELEMENT] = 0;				// The last element in the splits are removed
		SPLIT_ELEMENT -= 1;						// The number of splits is decreased by one

	};
	
	// ------------------------------------------------------------------------------------------------------------------------
    point_copy(G[0], SPLITTING_POINTS[0][0]);	// order-l point determined by T_{-}
    point_copy(G[1], SPLITTING_POINTS[1][0]);	// order-l point determined by T_{+}

    ec = lookup(L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1], key);	// To get current e_i in constant-time
    fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
    fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
    
	if( (isinfinity(G[0]) != 1) && (isequal(counter[L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1]], B[L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1]]) != 1) )
	{
		bc = isequal(ec >> 1, 0) & 1;					// Bit that determines if a dummy operation will be perfomed
		key[L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1]] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
        counter[L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1]] += 1;

		// Finally, we proceed by constructing the last one degree-(l_{NUMBER_OF_PRIMES[j] - 1}) isogeny (without evaluations)
		yISOG(tmp_Ts, current_A, SPLITTING_POINTS[(ec & 0x1)][0], current_A, L_STRATEGY[j][NUMBER_OF_PRIMES[j] - 1]);
	};

	// ------------------------------------------------------------------------------------------------------------------------
	point_copy(C, current_A);
};

/* ----------------------------------------------------------------------------------------------- *
   action_evaluation()
   inputs: a the secret key, the projective Edwards y-coordinates of y(T_{+}) = YT_{+}/ZT_{+}, the 
           Edwards curve constants A[0]:=a, and A[1]:=(a - d);
   output: the isogenous Edwards curve constants B[0]:=a' and B[1]:=(a' - d') determined by the action 
           evaluated at the secret key and public curve A
   NOTE: This is the algorithm that uses dummy operations
 * ----------------------------------------------------------------------------------------------- */
void action_evaluation(proj C, const uint8_t key[], const proj A)
{
	uint8_t tmp_key[N], counter[N];
	memcpy(tmp_key, key, sizeof(uint8_t) * N);

	proj current_A, current_Tp, current_Tm;
	point_copy(current_A, A);		// current Edwards curve constants a and (a -d)

	// Before constructing isogenies, we must to search for suitable point
	elligator(current_Tp, current_Tm, current_A);

	int current_i, current_j, local_j, local_i;
    uint8_t finished[N];                // flag that determines if the maximum number of isogeny constructions has been reached
    for(current_i = 0; current_i < N; current_i++)
        counter[current_i] = 0;

	for(current_i = 0; current_i < NUMBER_OF_DIFFERENT_STRATEGIES; current_i++)
	{
		for(current_j = 0; current_j < ROUNDS[current_i]; current_j++)
		{
			yDBL(current_Tp, current_Tp, current_A); yDBL(current_Tp, current_Tp, current_A);	// multiplication by [4]
			yDBL(current_Tm, current_Tm, current_A); yDBL(current_Tm, current_Tm, current_A);	// multiplication by [4]

			// Now, we multiply by l_j's until T'_{+} and T'_{-} have order at most the product of the first N[current_i] primes l_i's
			for(local_j = 0; local_j < (N - NUMBER_OF_PRIMES[current_i]) ; local_j++)
			{
				// Split corresponding to T_{+}
				yMUL(current_Tp, current_Tp, current_A, W_STRATEGY[current_i][local_j]);
				yMUL(current_Tm, current_Tm, current_A, W_STRATEGY[current_i][local_j]);
			};

			// Next, we make one single triangle evaluation
			strategy_evaluation(current_A, counter, tmp_key, current_Tp, current_Tm, current_A, current_i);

			// At this point, we must find another pair of points T_{+} and T_{-}. Thus, we use the projective elligator
			elligator(current_Tp, current_Tm, current_A);
		};
	};

	proj G[2], Z, K[(LARGE_L >> 1) + 1];
	int8_t ec = 0, bi = 0, ci = 0;
    uint32_t bc;
    while(1)
    {
        // number of degree-l_i isogenies don't reached
        local_j = 0;
        for(current_i = 0; current_i < N; current_i++)
        {
            ci = lookup(current_i, counter);
            bi = lookup(current_i, B);
            if( issmaller(ci, bi) == 0 )
                local_j += 1;
        };
        
        if(local_j == N)
            break;

        // Next, it is required to multiply the point by 4 and each l_i that doesn't belong to the current batch
        // T_{+}
        yDBL(current_Tp, current_Tp, current_A); // mult. by [2]
        yDBL(current_Tp, current_Tp, current_A); // mult. by [2]
        // T_{-}
        yDBL(current_Tm, current_Tm, current_A); // mult. by [2]
        yDBL(current_Tm, current_Tm, current_A); // mult. by [2]
        
        // number of degree-l_i isogenies don't reached
        local_j = 0;
        for(current_i = 0; current_i < N; current_i++)
        {
            ci = lookup(current_i, counter);
            bi = lookup(current_i, B);
            if( issmaller(ci, bi) == 0 )
            {
                yMUL(current_Tp, current_Tp, current_A, current_i); // Corresponding with T_{+}
                yMUL(current_Tm, current_Tm, current_A, current_i); // Corresponding with T_{-}
                
                finished[current_i] = 1;
            }
            else
            {
                finished[current_i] = 0;
                local_j += 1;
            }
        };
        
        local_i = 0;
        for(current_i = 0; current_i < N; current_i++)
        {
            if( finished[current_i] == 1 )
            { 
                //depends only on randomness
                continue;
            }
            else
            {
                // Now, a degree-(l_{batches[m][i]}) will be constructed. Let l = l_{batches[m][i]}.
                point_copy(G[0], current_Tm);	// order-l point determined by T_{-}
                point_copy(G[1], current_Tp);	// order-l point determined by T_{+}
                
                ec = lookup(current_i, tmp_key);			// To get current e_i in constant-time
				fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

                fp_cswap(current_Tm[0], current_Tp[0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                fp_cswap(current_Tm[1], current_Tp[1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                    
                for (current_j = (current_i + 1); current_j < N; current_j++)
                {
                    if( finished[current_j] == 0 )
                    {
                        //depends only on randomness
                        yMUL(G[0], G[0], current_A, current_j);	// Corresponding with T_{-}
                    };
                };
                
                if( local_i < (local_j - 1) )
                    yMUL(current_Tp, current_Tp, current_A, current_i);	// [l]T[1]

                if (isinfinity(G[0]) != 1)	// Depending on randomness
                {
                    bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed
                    
                    yISOG(K, current_A, G[0], current_A, current_i);
                    
                    if( local_i < (local_j - 1) )
                    {
                        yEVAL(current_Tm, current_Tm, K, current_i);	// evaluation of T[0]
                        yEVAL(current_Tp, current_Tp, K, current_i);	// evaluation of T[0]
                    }
                                        
                    tmp_key[current_i] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
                    counter[current_i] += 1;
                }

                fp_cswap(current_Tm[0], current_Tp[0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                fp_cswap(current_Tm[1], current_Tp[1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

                local_i += 1;
            };
        };
        

        // New Random torsion points T_{+} and T_{-}
        elligator(current_Tp, current_Tm, current_A);
    };
	point_copy(C, current_A);
};

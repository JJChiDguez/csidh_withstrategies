#include <math.h>
#include<time.h>

#include "fp.h"
#include "edwards_curve.h"

// A utility function to swap two elements
void swap(uint64_t* a, uint64_t* b)
{
	uint64_t t = *a;
	*a = *b;
	*b = t;
}

/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition (uint64_t arr[], int low, int high)
{
	uint64_t pivot = arr[high];    // pivot
	int i = (low - 1);  // Index of smaller element

	int j;
	for(j = low; j <= high- 1; j++)
	{
		// If current element is smaller than the pivot
		if(arr[j] < pivot)
		{
			i++;    // increment index of smaller element
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}

/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
void quicksort(uint64_t arr[], int low, int high)
{
    	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now at right place */
		int pi = partition(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quicksort(arr, low, pi - 1);
		quicksort(arr, pi + 1, high);
	}
}

// Measuring the perfomance
static uint64_t get_cycles()
{
	uint32_t lo, hi;
   	asm volatile("rdtsc":"=a"(lo),"=d"(hi));
   	return ((uint64_t)hi<<32) | lo;
};

static uint8_t csidh(proj out, const uint8_t sk[], const proj in)
{
	if (!validate(in)) {
		return 0;
	};
	action_evaluation(out, sk, in);
	return 1;
};

unsigned long its = 1024;

int main()
{
	unsigned int i;

	uint64_t add_min = 0xFFFFFFFFFFFFFFFF, add_max = 0, 
	         sqr_min = 0xFFFFFFFFFFFFFFFF, sqr_max = 0, 
	         mul_min = 0xFFFFFFFFFFFFFFFF, mul_max = 0;

	float add_median = 0, add_mean = 0, add_variance = 0,
	      sqr_median = 0, sqr_mean = 0, sqr_variance = 0,
	      mul_median = 0, mul_mean = 0, mul_variance = 0;

	uint64_t add_sample[its],
	         sqr_sample[its],
	         mul_sample[its];

	// ---
	uint8_t key[N];
	proj random_E;
	point_copy(random_E, E);

	for(i = 0; i < its; ++i)
	{

		if (its >= 100 && i % (its / 100) == 0) {
			printf("Doing %lu iterations of action with validation key:\t", its);
            printf("%2lu%%", 100 * i / its);
            fflush(stdout);
            printf("\r\x1b[K");
        }

		random_key(key);
		init_counters();        // counters of additions, squarings, and multiplications are set as zero
		assert(csidh(random_E, key, random_E));
		
		// ---
		add_sample[i] = fpadd;
		sqr_sample[i] = fpsqr;
		mul_sample[i] = fpmul;

		/**************************************/
		if(add_min > add_sample[i])
			add_min = add_sample[i];
		if(sqr_min > sqr_sample[i])
			sqr_min = sqr_sample[i];
		if(mul_min > mul_sample[i])
			mul_min = mul_sample[i];

		/**************************************/
		if(add_max < add_sample[i])
			add_max = add_sample[i];
		if(sqr_max < sqr_sample[i])
			sqr_max = sqr_sample[i];
		if(mul_max < mul_sample[i])
			mul_max = mul_sample[i];
		
		/**************************************/
		add_mean += (float)add_sample[i];
		sqr_mean += (float)sqr_sample[i];
		mul_mean += (float)mul_sample[i];
	};


	add_mean = add_mean / ((float)its * 1.0);
	sqr_mean = sqr_mean / ((float)its * 1.0);
	mul_mean = mul_mean / ((float)its * 1.0);

	for (i = 0; i < its; ++i)
	{
		add_variance += (add_sample[i] - add_mean)*(add_sample[i] - add_mean);
		sqr_variance += (sqr_sample[i] - sqr_mean)*(sqr_sample[i] - sqr_mean);
		mul_variance += (mul_sample[i] - mul_mean)*(mul_sample[i] - mul_mean);
	};

	add_variance = add_variance / ((float)its - 1.0);
	sqr_variance = sqr_variance / ((float)its - 1.0);
	mul_variance = mul_variance / ((float)its - 1.0);
	
	quicksort(mul_sample, 0, its - 1);
	quicksort(sqr_sample, 0, its - 1);
	quicksort(add_sample, 0, its - 1);

	add_median = (float)(add_sample[its/2] + add_sample[its/2 - 1]) / 2.0;
	sqr_median = (float)(sqr_sample[its/2] + sqr_sample[its/2 - 1]) / 2.0;
	mul_median = (float)(mul_sample[its/2] + mul_sample[its/2 - 1]) / 2.0;

	printf("\x1b[01;33mIterations: %lu\x1b[0m\n\n", its);

	printf("\x1b[33mMedian costs:\x1b[0m\n");
	printf("\t %f additions,\n", add_median);
	printf("\t\x1b[32m %f squarings,\x1b[0m\n", sqr_median);
	printf("\t\x1b[31m %f multiplications.\x1b[0m\n", mul_median);

	printf("\n");

	printf("\x1b[33mAverage costs:\x1b[0m\n");
	printf("\t %f additions,\n", add_mean);
	printf("\t\x1b[32m %f squarings,\x1b[0m\n", sqr_mean);
	printf("\t\x1b[31m %f multiplications.\x1b[0m\n", mul_mean);

	printf("\n");

	printf("\x1b[33mStandard deviation of the costs:\x1b[0m\n");
	printf("\t %f additions,\n", sqrt(add_variance));
	printf("\t\x1b[32m %f squarings,\x1b[0m\n", sqrt(sqr_variance));
	printf("\t\x1b[31m %f multiplications.\x1b[0m\n", sqrt(mul_variance));

	printf("\n");

	printf("\x1b[33mMinimum costs:\x1b[0m\n");
	printf("\t %lu additions,\n", add_min);
	printf("\t\x1b[32m %lu squarings,\x1b[0m\n", sqr_min);
	printf("\t\x1b[31m %lu multiplications.\x1b[0m\n", mul_min);

	printf("\n");

	printf("\x1b[33mMaximum costs:\x1b[0m\n");
	printf("\t %lu additions,\n", add_max);
	printf("\t\x1b[32m %lu squarings,\x1b[0m\n", sqr_max);
	printf("\t\x1b[31m %lu multiplications.\x1b[0m\n", mul_max);

	return 0;
};

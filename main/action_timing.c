#include <math.h>
#include<time.h>

#include "fp.h"
#include "edwards_curve.h"

// A utility function to swap two elements
void swap(double* a, double* b)
{
        double t = *a;
        *a = *b;
        *b = t;
}

/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition (double arr[], int low, int high)
{
        double pivot = arr[high];    // pivot
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
void quicksort(double arr[], int low, int high)
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

	double cc_min = 0xFFFFFFFFFFFFFFFF, cc_max = 0;
	double cc_median = 0, cc_mean = 0, cc_variance = 0;

	double cc_sample[its], cc_0, cc_1;

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
		cc_0 = get_cycles();
		assert(csidh(random_E, key, random_E));
		cc_1 = get_cycles();
		
		cc_sample[i] = (double)(cc_1 - cc_0) / (1000000.0) ;
		// ---


		/**************************************/
		if(cc_min > cc_sample[i])
			cc_min = cc_sample[i];

		/**************************************/
		if(cc_max < cc_sample[i])
			cc_max = cc_sample[i];
		
		/**************************************/
		cc_mean += (double)cc_sample[i];
	};


	cc_mean = cc_mean / ((double)its * 1.0);

	for (i = 0; i < its; ++i)
	{
		cc_variance += (cc_sample[i] - cc_mean)*(cc_sample[i] - cc_mean);
	};

	cc_variance = cc_variance / ((double)its - 1.0);
	
	quicksort(cc_sample, 0, its - 1);
	cc_median = (cc_sample[its/2] + cc_sample[its/2 - 1]) / 2.0;

	printf("\x1b[01;33mIterations: %lu\x1b[0m\n\n", its);

	printf("\x1b[33mMedian number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_median);
	printf("\x1b[33mAverage number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_mean);
	printf("\x1b[33mStandar deviation of the number of clock cycles: \x1b[32m %f \x1b[0m\n", sqrt(cc_variance));
	printf("\x1b[33mMinimum of the number clock cycles: \x1b[32m %f \x1b[0m\n", cc_min);
	printf("\x1b[33mMaximum of the number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_max);
	printf("\n");

	return 0;
};

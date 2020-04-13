# This is Python3-code implementation of CSIDH protocol.

## Execution

### CSIDH protocol and Benchmark

- python3 csidh.py p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY]
- python3 bench.py p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY] NUMBER_OF_RUNS


Here, p[NUMBER_OF_BITS] must contain the list of Small Odd Primes: l_1 l_2 l_3 ... l_{n - 2} l_{n - 1} l_n
However, in this case only the 511-bit prime of CSIDH-512 has been used.

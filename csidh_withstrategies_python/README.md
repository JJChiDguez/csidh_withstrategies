# This is Python3-code implementation of CSIDH protocol.

## Execution

### CSIDH protocol and Benchmark

- python3 csidh.py p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY]
- python3 bench.py p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY] NUMBER_OF_RUNS

Here, p[NUMBER_OF_BITS] must contain the list of Small Odd Primes: l_1 l_2 l_3 ... l_{n - 2} l_{n - 1} l_n

### Example of a run

However, experiments were ran with the 511-bit and 1024-bit primes used in CSIDH-512 CSIDH-1024 protocol implementations given in https://csidh.isogeny.org/index.html.

- python3 csidh.py p512 wd2 STRATEGY
- python3 bench.py p512 wd2 STRATEGY 1024

- python3 csidh.py p1024 wd2 STRATEGY
- python3 bench.py p1024 wd2 STRATEGY 1024

# Python3-code implementation of CSIDH protocol.

This implementation corresponds with the work titled _"Optimal strategies for CSIDH"_ by **Jesús-Javier Chi-Domínguez** and **Francisco Rodríguez-Henríquez**.

## Dependencies

The provided Python3-code requires the package `python3-progress`, `python3-numpy`. Thus, these package can be installed (in Ubuntu distro) by running the next commands

```bash
sudo apt-get install python3-progress
sudo apt-get install python3-numpy
```

## Execution

### CSIDH protocol, Benchmark, and suitable bounds search

- python3 csidh.py p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY]
- python3 bench.py p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY] NUMBER_OF_RUNS
- python3 suitable_bounds.py p[NUMBER_OF_BITS] [wd1/wd2/df] STRATEGY

Here, p[NUMBER_OF_BITS] must contain the list of Small Odd Primes: l_1 l_2 l_3 ... l_{n - 2} l_{n - 1} l_n

### Example of a run

However, experiments were ran with the 511-bit and 1024-bit primes used in CSIDH-512 CSIDH-1024 protocol implementations given in https://csidh.isogeny.org/index.html.

```bash
# CSIDH-512
python3 csidh.py p512 wd2 STRATEGY
python3 bench.py p512 wd2 STRATEGY 1024
# CSIDH-1024
python3 csidh.py p1024 wd2 STRATEGY
python3 bench.py p1024 wd2 STRATEGY 1024

# Looking for suitable bounds (CSIDH-512)
python3 suitable_bounds.py p512 wd2 STRATEGY
python3 suitable_bounds.py p512 wd1 STRATEGY
python3 suitable_bounds.py p512 df STRATEGY
# Looking for suitable bounds (CSIDH-1024)
python3 suitable_bounds.py p1024 wd2 STRATEGY
python3 suitable_bounds.py p1024 wd1 STRATEGY
python3 suitable_bounds.py p1024 df STRATEGY
```

### Remarks
- The integer bound vector required in the files `csidh.py` and `bench.py` must be in ascending way (going from small to large values of `l_i`); this is because, our strategy method computes degree-`l_{n - i}` isogenies at each step.
- The vector output given by the file `suitable_bounds.py` are presented in descending way (going from large to small values of `l_i`)
- The file `suitable_bounds.py` gives an expected cost without taking in count the cost of the _missing part_ (unreached `l_i's`)
- Recall, `suitable_bounds.py` implements a greedy algorithm, and its time execution depends on the initial bound value `m` 

## Funding

This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 804476).
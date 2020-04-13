from progress.bar import Bar
from functools import reduce
import sys
import random
from math import ceil, floor, log
import numpy as np
import statistics


bitlength = lambda x: len(bin(x)[2:])                       # number of bits
hamming_weight = lambda x: bin(x).count("1");               # hamming weight: number of bits equal 1

sign = lambda x: (1, -1)[x < 0]                             # Sign of an integer
isequal = { True : 1 , False : 0 }                          # Simulating constant-time integer comparison

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# Inputs

if( ((len(sys.argv) != 5) and (sys.argv[0] == 'bench.py')) or ((len(sys.argv) != 4) and (sys.argv[0] == 'csidh.py')) ):
    print("[ERROR]\tThis script receives two inputs:")
    if (sys.argv[0] == 'bench.py'):
        print("\t-] python3 %s p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY] NUMBER_OF_RUNS" %sys.argv[0])
    else:
        print("\t-] python3 %s p[NUMBER_OF_BITS] [wd1/wd2/df] [SIMBA/SIMBA_WITH_STRATEGY/STRATEGY]" %sys.argv[0])
    print("Here, p[NUMBER_OF_BITS] must contain the list of Small Odd Primes: l_1 l_2 l_3 ... l_{n - 2} l_{n - 1} l_n")
    exit(1)

try:

    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
    f = open('./sop/' + sys.argv[1])
    L = f.read()
    L = [ int(l) for l in L.split() ]
    n = len(L)
    f.close()

except IOError:

    print("[ERROR]\tFile not accessible")
    exit(2)

TYPE = sys.argv[2]                                      # wd1, wd2, or dummy-free
APPROACH = sys.argv[3]                                  # SIMBA, SIMBA with strategies, or Strategies

if( (APPROACH != 'STRATEGY') and (APPROACH != 'SIMBA') and (APPROACH != 'SIMBA_WITH_STRATEGY') ):
    print("[ERROR]\tApproach must be one of the followings:")
    print("\t-] SIMBA")
    print("\t-] SIMBA_WITH_STRATEGY")
    print("\t-] STRATEGY")
    exit()

if(sys.argv[0] == 'bench.py'):
    MAX_iter = int(sys.argv[4])

p = 4 * reduce(lambda x,y : (x*y), L) - 1               # p := 4 * l_0 * ... * l_n - 1
p_minus_one_halves = (p - 1) // 2                       # (p - 1) / 2

bits_of_4sqrt_of_p = ceil(bitlength(p) / 2.0) + 2
# --------------------------------------------------------------------------------------------------------------------------------
# Checking if p is composite

def _try_composite(a, d, n, s):
    if pow(a, d, n) == 1:
        return False
    for i in range(s):
        if pow(a, 2**i * d, n) == n-1:
            return False
    return True # n  is definitely composite

def is_prime(n):
    """
    Miller-Rabin primality test.
 
    A return value of False means n is certainly not prime. A return value of
    True means n is very likely a prime.
    """
    if n!=int(n):
        return False
    n=int(n)
    #Miller-Rabin test for prime
    if n==0 or n==1 or n==4 or n==6 or n==8 or n==9:
        return False
 
    if n==2 or n==3 or n==5 or n==7:
        return True
    s = 0
    d = n-1
    while d%2==0:
        d>>=1
        s+=1
    assert(2**s * d == n-1)
 
    def trial_composite(a):
        if pow(a, d, n) == 1:
            return False
        for i in range(s):
            if pow(a, 2**i * d, n) == n-1:
                return False
        return True  
 
    for i in range(128):	#number of trials 
        a = random.randrange(2, n)
        if trial_composite(a):
            return False
 
    return True

if is_prime(p) == False:
	print("[ERROR]\tThe integer number p := 4 * l_1, * ... * l_n - 1 is not prime where L := ", L)
	exit(3)

if( (sys.argv[0] != 'header.py') ):
    print("/*")
    print("The prime number to be used is\n")
    print("           %03d" % n)
    print("        ~~~~~~~~~")
    print("         |     | ")
    print("p := 4 x |     | l_j   -   1, where each l_j is a small odd prime ")
    print("           j=1  ")
    print("*/\n")

# At this point, p is very likely a prime. Thus, we can continue
# --------------------------------------------------------------------------------------------------------------------------------

# Jacobi symbol used for checking if an integer has square-root in fp
def jacobi(a, n):

    assert(n > a > 0 and n%2 == 1)
    t = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t
        a, n = n, a
        if a % 4 == n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    else:
        return 0

# Extended GCD
def xgcd(aa, bb):
	lastremainder, remainder = abs(aa), abs(bb)
	x, lastx, y, lasty = 0, 1, 1, 0
	while remainder:
		lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
		x, lastx = lastx - quotient*x, x
		y, lasty = lasty - quotient*y, y

	return lastremainder, lastx * (-1 if aa < 0 else 1), lasty * (-1 if bb < 0 else 1)

# Modular inverse
def fp_inv(a):
	g, x, y = xgcd(a, p)
	#if g != 1:
	#	raise ValueError
	return x % p

# Modular addition
def fp_add(a, b):
    return (a + b) % p

# Modular substraction
def fp_sub(a, b):
    return (a - b) % p

# Modular multiplication
def fp_mul(a, b):
    return (a * b) % p

# Modular squaring
def fp_sqr(a):
    return (a ** 2) % p

# constant-time swap
def fp_cswap(x, y, b):
    z = list([x, y])
    z = list(z[::(1 - 2*b)])
    return z[0], z[1]

# Modular exponentiation
def fp_exp(a, e):
    b = 1
    while e > 0:
        if (e & 0) == 1:
            b = fp_mul(a, b)
        a = fp_sqr(a)
        e = e >> 1

    return b

# --------------------------------------------------------------------------------------------------------------------------------
# The next variable is used for counting the number of field operations
FIELD_OPERATIONS_PERFORMED = np.array([0.0, 0.0, 0.0])    # (M, S, a)

def SET_FIELD_OPERATIONS_TO_ZERO():

    # The field operations are set to 0 (subtracting by itself)
    global FIELD_OPERATIONS_PERFORMED
    FIELD_OPERATIONS_PERFORMED -= FIELD_OPERATIONS_PERFORMED

'''
    chunks()
    inputs: a string, a list, and the maximum  number of elements in each chunk
    -----
    NOTE: This function divide the input list into len(L) / k chunks.
'''
chunks = lambda NAME, L, n : [NAME + ' =\t{'] +\
                             [ '\t' + ', '.join(list(map(format, L[i * n:(i + 1) * n], ['3d']*n))) for i in range((len(L) + n - 1) // n )] +\
                             ['\t};']
'''
    printl()
    inputs: a string, a list, and the maximum number k of elements in each chunk
    -----
    NOTE: this function prints a given list by chunks of size k.
'''
def printl(NAME, L, k):

    to_print = chunks(NAME, L, k)
    print(to_print[0])
    for i in range(1, len(to_print) - 2):
        print(to_print[i] + ",")

    print(to_print[len(to_print) - 2])
    print(to_print[len(to_print) - 1])

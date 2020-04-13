from fp import *
from montgomery import *
#from isogeny import *
from isogeny_hybrid import *

if TYPE == 'wd1':
    # Using one torsion point and dummy isogeny constructions
    from gae_wd1 import *
elif TYPE == 'wd2':
    # Using two torsion points and dummy isogeny constructions
    from gae_wd2 import *
elif TYPE == 'df':
    # Dummy-free approach by using two torsion points
    from gae_df import *
else:
    print("[ERROR]\tInvalid TYPE. Second input must be one of the following:")
    print("\t-] wd1 <- using ONE torsion point and dummy isogeny constructions,")
    print("\t-] wd2 <- using TWO torsion points and dummy isogeny constructions, and")
    print("\t-] df <- using a dummy-free approach and TWO torsion points.")
    exit(7)

''' -------------------------------------------------------------------------------------
    Framework
    ------------------------------------------------------------------------------------- '''
    
print("\n")
print("p := 0x%X;" % p)
print("fp := GF(p);")
print("P<x> := PolynomialRing(fp);");
print("fp2<i> := ext<fp | x^2 + 1>;")
print("P<x> := PolynomialRing(fp2);");

A = [2, 4]
print("E_A := EllipticCurve(x^3 + 0x%X * x^2 + x);\n" % coeff(A))


''' -------------------------------------------------------------------------------------
    Number of degree-(l_i) isogeny constructions to be performed: m_i
    \sigma and \kappa parameters of SIMBA method
    ------------------------------------------------------------------------------------- '''

# ==========================================================================

if( sys.argv[1] == 'p512' ):

    if(TYPE == 'wd1'):
        # ====== [MCR style, CSIDH-512] each m_i corresponds with the given in MCR18
        #m = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5]

        # ===== Suitable bounds for this work
        m = [15, 18, 20, 21, 21, 22, 22, 22, 22, 22, 22, 19, 20, 22, 23, 23, 23, 23, 23, 23, 23, 21, 23, 20, 16, 16, 16, 15, 14, 12, 13, 12, 11, 11, 10, 10, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3]

        # ====== [MCR style, CSIDH-512] case when each m_i is equal to 10, and it implies a key space of (10 + 1)^74 ~ 2^256 ~ p^1/4
        #m = [10] * n
        
        #sigma, kappa = 1, 10 # when only one strategy is required; that is, m = (10, 10, ..., 10)
        sigma, kappa = 5, 11 # MCR & dummy-free [The given one in MCR18] CSIDH-512
    
    if(TYPE == 'wd2'):
        # ====== [OAYT style, CSIDH-512] each m_i corresponds with the given in OAYT19
        #m = [5, 6, 7, 7, 7, 7, 7, 8, 8, 8, 9, 10, 10, 10, 10, 9, 9, 9, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1]

        # ===== Suitable bounds for this work
        m = [7, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 9, 11, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]

        # ====== [OAYT style, CSIDH-512] case when each m_i is equal to 5, and it implies a key space of (2*5 + 1)^74 ~ 2^256 ~ p^1/4
        #m = [5] * n
        
        #sigma, kappa = 1, 5 # when only one strategy is required; that is, m = (5, 5, ..., 5)
        sigma, kappa = 3, 8 # OAYT [The given one in OAYT19] CSIDH-512
    
    if(TYPE == 'df'):
        # ====== [dummy-free style, CSIDH-512] each m_i corresponds with the given in MCR18 (it is the same as MCR style)
        #m = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5]

        # ===== Suitable bounds for this work
        m = [15, 18, 20, 21, 21, 22, 22, 22, 22, 22, 22, 19, 20, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 19, 16, 16, 16, 15, 14, 12, 13, 12, 11, 11, 11, 9, 9, 9, 9, 8, 8, 8, 8, 7, 8, 6, 6, 6, 6, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3]

        # ====== [dummy-free style, CSIDH-512] case when each m_i is equal to 10, and it implies a key space of (10 + 1)^74 ~ 2^256 ~ p^1/4
        #m = [10] * n
        
        #sigma, kappa = 1, 10 # when only one strategy is required; that is, m = (10, 10, ..., 10)
        sigma, kappa = 5, 11 # MCR & dummy-free [The given one in MCR18] CSIDH-512
        
else:
    print("[ERROR]\tMissing bound. To set the maximum number of isogeny constructions, and sigma and kappa from SIMBA method (to add they in ELSE statement in line 85 of file csidh.py)")
    exit(11)

# ==========================================================================

if len(set(m)) > 1:
    # Maximum number of degree-(l_i) isogeny constructions is m_i (different for each l_i)
    LABEL_m = ''.join(['{0:0{1}d}'.format(m_i,2) for m_i in m ])
else:
    # Maximum number of degree-(l_i) isogeny constructions is m (the same for each l_i)
    LABEL_m = ''.join(['{0:0{1}d}'.format(m_i,2) for m_i in set(m) ])

try:

    # List of Small Odd Primes, L := [l_0, ..., l_{n-1}]
    m_prime = [ geometric_serie(m[k], L[k]) for k in range(n) ]
    r_out, L_out, R_out = rounds(m_prime[::-1], n)
    for j in range(0, len(r_out), 1):

        R_out[j] = list([L[::-1][k] for k in R_out[j]])
        L_out[j] = list([L[::-1][k] for k in L_out[j]])

    f = open('./strategies/' + sys.argv[1] + '_' + TYPE + '_' + str(SQR) + '_' + str(ADD) + '_' + LABEL_m)
    print("\"Strategies to be read from a file\";\n")
    S_out = []
    for i in range(0, len(r_out), 1):

        tmp = f.readline()
        tmp = [ int(b) for b in tmp.split() ]
        S_out.append(tmp)

    f.close()

except IOError:

    print("\"Strategies to be computed\";\n")
    C_out, L_out, R_out, S_out, r_out = strategy_block_cost(L[::-1], m[::-1])
    f = open('./strategies/' + sys.argv[1] + '_' + TYPE + '_' + str(SQR) + '_' + str(ADD) + '_' + LABEL_m,'w')
    for i in range(0, len(r_out)):

        f.writelines(' '.join([ str(tmp) for tmp in S_out[i]]) + '\n')

    f.close()

print("// All the experiments are assuming S = %1.6f x M and a = %1.6f x M. The measures are given in millions of field operations.\n" % (SQR, ADD))

''' -------------------------------------------------------------------------------------
    Main
    ------------------------------------------------------------------------------------- '''

print("// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
print("// Public Key Generation")
print("// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

temporal_m = list(set(m))
T_p, T_m = full_torsion_points(A)

# ------------------------------------------------------------------------- Alice
SET_FIELD_OPERATIONS_TO_ZERO()
public_validation = validate(A)
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)
print("// Number of field operations (validate):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

a_private = random_key(m)
SET_FIELD_OPERATIONS_TO_ZERO()
if('SIMBA' in APPROACH):
    a_public = SIMBA(A, a_private, L, n, sigma, kappa, m, APPROACH == 'SIMBA_WITH_STRATEGY')
else:
    if len(temporal_m) == 1:
        if temporal_m[0] > 1:
            a_public = GAE(A, a_private, [L_out[0]], [ [] ], [S_out[0]], temporal_m, m)
        else:
            if TYPE == 'wd1':
                # ONE torsion  point required
                a_public, m_tmp, a_tmp = evaluate_strategy(A, T_p, L_out[0], S_out[0], n, m, a_private)
            else:
                # TWO torsion points required
                a_public, m_tmp, a_tmp = evaluate_strategy(A, list([list(T_m), list(T_p)]), L_out[0], S_out[0], n, m, a_private)
    else:
        a_public = GAE(A, a_private, L_out, R_out, S_out, r_out, m)
    
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)

print("// Public key corresponding to Alice")
print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
print("a_public := EllipticCurve(x^3 + 0x%X * x^2 + x);\n" % coeff(a_public))

# ------------------------------------------------------------------------- Bob
SET_FIELD_OPERATIONS_TO_ZERO()
public_validation = validate(A)
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)
print("// Number of field operations (validate):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

b_private = random_key(m)
SET_FIELD_OPERATIONS_TO_ZERO()
if('SIMBA' in APPROACH):
    b_public = SIMBA(A, b_private, L, n, sigma, kappa, m, APPROACH == 'SIMBA_WITH_STRATEGY')
else:
    if len(temporal_m) == 1:
        if temporal_m[0] > 1:
            b_public = GAE(A, b_private, [L_out[0]], [ [] ], [S_out[0]], temporal_m, m)
        else:
            if TYPE == 'wd1':
                # ONE torsion  point required
                b_public, m_tmp, b_tmp = evaluate_strategy(A, T_p, L_out[0], S_out[0], n, m, b_private)
            else:
                # TWO torsion points required
                b_public, m_tmp, b_tmp = evaluate_strategy(A, list([list(T_m), list(T_p)]), L_out[0], S_out[0], n, m, b_private)
    else:
        b_public = GAE(A, b_private, L_out, R_out, S_out, r_out, m)
    
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)

print("// Public key corresponding to Bob")
print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
print("b_public := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(b_public))

print("\n// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
print("// Secret Sharing Computation")
print("// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

# ------------------------------------------------------------------------- Alice
SET_FIELD_OPERATIONS_TO_ZERO()
public_validation = validate(b_public)
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)
print("// Number of field operations (validate):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

Tp_b, Tm_b = full_torsion_points(b_public)
SET_FIELD_OPERATIONS_TO_ZERO()
if('SIMBA' in APPROACH):
    ss_a = SIMBA(b_public, a_private, L, n, sigma, kappa, m, APPROACH == 'SIMBA_WITH_STRATEGY')
else:
    if len(temporal_m) == 1:
        if temporal_m[0] > 1:
            ss_a = GAE(b_public, a_private, [L_out[0]], [ [] ], [S_out[0]], temporal_m, m)
        else:
            if TYPE == 'wd1':
                # ONE torsion  point required
                ss_a, m_tmp, a_tmp = evaluate_strategy(b_public, Tp_b, L_out[0], S_out[0], n, m, a_private)
            else:
                # TWO torsion points required
                ss_a, m_tmp, a_tmp = evaluate_strategy(b_public, list([list(Tm_b), list(Tp_b)]), L_out[0], S_out[0], n, m, a_private)
    else:
        ss_a = GAE(b_public, a_private, L_out, R_out, S_out, r_out, m)
    
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)

print("// Public key corresponding to Alice")
print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
print("ss_a := EllipticCurve(x^3 + 0x%X * x^2 + x);\n" % coeff(ss_a))

# ------------------------------------------------------------------------- Bob
SET_FIELD_OPERATIONS_TO_ZERO()
public_validation = validate(a_public)
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)
print("// Number of field operations (validate):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )

Tp_a, Tm_a = full_torsion_points(a_public)
SET_FIELD_OPERATIONS_TO_ZERO()
if('SIMBA' in APPROACH):
    ss_b = SIMBA(a_public, b_private, L, n, sigma, kappa, m, APPROACH == 'SIMBA_WITH_STRATEGY')
else:
    if len(temporal_m) == 1:
        if temporal_m[0] > 1:
            ss_b = GAE(a_public, b_private, [L_out[0]], [ [] ], [S_out[0]], temporal_m, m)
        else:
            if TYPE == 'wd1':
                # ONE torsion  point required
                ss_b, m_tmp, b_tmp = evaluate_strategy(a_public, Tp_a, L_out[0], S_out[0], n, m, b_private)
            else:
                # TWO torsion points required
                ss_b, m_tmp, b_tmp = evaluate_strategy(a_public, list([list(Tm_a), list(Tp_a)]), L_out[0], S_out[0], n, m, b_private)
    else:
        ss_b = GAE(a_public, b_private, L_out, R_out, S_out, r_out, m)
    
RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)

print("// Public key corresponding to Bob")
print("// Number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
print("ss_b := EllipticCurve(x^3 + 0x%X * x^2 + x);" % coeff(ss_b))

print("// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
print("// Private Key Generation")
print("// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
print("// Private key corresponding to Alice")
print("/*")
printl("a_private", a_private, n // 7)
print("*/")
print("// Private key corresponding to Bob")
print("/*")
printl("b_private", b_private, n // 7)
print("*/")
print("// Maximum number of degree-(\ell_i) isogeny constructions: m_i")
print("/*")
printl("m", m[::-1], n // 7)
print("*/")
print("\n// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
if( coeff(ss_a) == coeff(ss_b) ):
    print('\x1b[0;30;43m' + '\"Successfully passed!\";' + '\x1b[0m')
else:
    print('\x1b[0;30;41m' + '\"Great Scott!... The sky is falling. NOT PASSED!!!\"' + '\x1b[0m')


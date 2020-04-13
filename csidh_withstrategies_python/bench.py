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
    print("[ERROR]\tMissing bound. To set the maximum number of isogeny constructions, and sigma and kappa from SIMBA method (to add they in ELSE statement in line 85 of file bench.py)")
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

''' -------------------------------------------------------------------------------------
    Main
    ------------------------------------------------------------------------------------- '''
print("// All the experiments are assuming S = %1.6f x M and a = %1.6f x M. The measures are given in millions of field operations.\n" % (SQR, ADD))

print("// Maximum number of degree-(\ell_i) isogeny constructions: m_i")
print("/*")
printl("m", m, n // 7)
print("*/")

if( APPROACH != "STRATEGY"):

    SAMPLE = [ [0.0, 0.0, 0.0] ] * MAX_iter
    SAMPLE_VALIDATE = [ [0.0, 0.0, 0.0] ] * MAX_iter
    bar = Bar('Experiments corresponding to SIMBA-' + str(sigma) + '-' + str(kappa) , max=MAX_iter)
    B = list(A)
    for main_i in range(MAX_iter):

        e = random_key(m)
        
        SET_FIELD_OPERATIONS_TO_ZERO()
        B = SIMBA(B, e, L, n, sigma, kappa, m, APPROACH == "SIMBA_WITH_STRATEGY")
        SAMPLE[main_i] = list(FIELD_OPERATIONS_PERFORMED)
        #print("Random(EllipticCurve(x^3 + 0x%X * x^2 + x)) * (p+1);" % coeff(b))

        SET_FIELD_OPERATIONS_TO_ZERO()
        V = validate(B)
        SAMPLE_VALIDATE[main_i] = list(FIELD_OPERATIONS_PERFORMED)
                
        bar.next()
    bar.finish()

    AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE ]), statistics.mean([ ops[1] for ops in SAMPLE ]), statistics.mean([ ops[2] for ops in SAMPLE ]) ]
    print("Average number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )

else:

    e = random_key(m)
    B = list(A)
    B = SIMBA(B, e, L, n, sigma, kappa, m, APPROACH == "SIMBA_WITH_STRATEGY")
    
    tmp = list(set(m))
    SAMPLE = [ [0.0, 0.0, 0.0] ] * MAX_iter
    SAMPLE_FULL_TORSION_POINT = [ [0.0, 0.0, 0.0] ] * MAX_iter
    SAMPLE_VALIDATE = [ [0.0, 0.0, 0.0] ] * MAX_iter
    bar = Bar('Experiments corresponding to this work' , max=MAX_iter)

    for main_i in range(MAX_iter):

        e = random_key(m)
        
        if len(tmp) == 1:
        
            if tmp[0] > 1:
                # NOT fully constant-time group action evaluation
                SET_FIELD_OPERATIONS_TO_ZERO()
                B = GAE(B, e, [L_out[0]], [ [] ], [S_out[0]], tmp, m)
                SAMPLE[main_i] = list(FIELD_OPERATIONS_PERFORMED)
            else:
                
                # Fully constant-time group action evaluation (only one strategy evaluation is required)
                SET_FIELD_OPERATIONS_TO_ZERO()
                T_p, T_m = full_torsion_points(B)
                SAMPLE_FULL_TORSION_POINT[main_i] = list(FIELD_OPERATIONS_PERFORMED)
            
        else:
            
            SET_FIELD_OPERATIONS_TO_ZERO()
            B = GAE(B, e, L_out, R_out, S_out, r_out, m)
            SAMPLE[main_i] = list(FIELD_OPERATIONS_PERFORMED)
            
        SET_FIELD_OPERATIONS_TO_ZERO()
        V = validate(B)
        SAMPLE_VALIDATE[main_i] = list(FIELD_OPERATIONS_PERFORMED)
        
        #print("Random(EllipticCurve(x^3 + 0x%X * x^2 + x)) * (p+1);" % coeff(B))
        
        bar.next()
    bar.finish()

    # Fully constant-time group action evaluation (only one strategy evaluation is required)
    if (len(tmp) == 1) and (tmp[0] == 1):

        SET_FIELD_OPERATIONS_TO_ZERO()
        #B, m_tmp, e_tmp = evaluate_strategy(B, T_p, L_out[0], S_out[0], n, m, e)
        if TYPE == 'wd1':
            # ONE torsion  point required
            B, m_tmp, e_tmp = evaluate_strategy(B, T_p, L_out[0], S_out[0], n, m, e)
        else:
            # TWO torsion points required
            B, m_tmp, e_tmp = evaluate_strategy(B, list([list(T_m), list(T_p)]), L_out[0], S_out[0], n, m, e)
            
        RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)
        print("Fixed number of field operations (GAE):\t\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
        
        AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE_FULL_TORSION_POINT ]), statistics.mean([ ops[1] for ops in SAMPLE_FULL_TORSION_POINT ]), statistics.mean([ ops[2] for ops in SAMPLE_FULL_TORSION_POINT ]) ]
        print("Randomly selection of full torsion points:\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )
    else:
        AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE ]), statistics.mean([ ops[1] for ops in SAMPLE ]), statistics.mean([ ops[2] for ops in SAMPLE ]) ]
        print("Average number of field operations (GAE):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )

AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE_VALIDATE ]), statistics.mean([ ops[1] for ops in SAMPLE_VALIDATE ]), statistics.mean([ ops[2] for ops in SAMPLE_VALIDATE ]) ]
print("Average number of field operations (validate):\t%1.6f x M + %1.6f x S + %1.6f x a := %1.6f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )

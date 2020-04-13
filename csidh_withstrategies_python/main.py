from fp import *
from montgomery import *
#from isogeny import *
from isogeny_hybrid import *

if TYPE == 'wd1':
    from gae_wd1 import *
elif TYPE == 'wd2':
    print("Not implemented yet!!!")
    exit(7)
elif TYPE == 'df':
    print("Not implemented yet!!!")
    exit(7)
else:
    print("[ERROR]\tInvalid TYPE. Second input must be one of the following:")
    print("\t-] wd1 <- using ONE torsion point and dummy isogeny constructions,")
    print("\t-] wd2 <- using TWO torsion points and dummy isogeny constructions, and")
    print("\t-] df <- using a dummy-free approach and TWO torsion points.")
    exit(7)

print("\n")
print("p := 0x%X;" % p)
print("fp := GF(p);")
print("P<x> := PolynomialRing(fp);");
print("fp2<i> := ext<fp | x^2 + 1>;")
print("P<x> := PolynomialRing(fp2);");

A = [2, 4]
print("E_A := EllipticCurve(x^3 + %d * x^2 + x);\n" % coeff(A))

#m = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5]
m = [1] * n

if len(set(m)) > 1:
    LABEL_m = ''.join(['{0:0{1}d}'.format(m_i,2) for m_i in m ])
else:
    LABEL_m = ''.join(['{0:0{1}d}'.format(m_i,2) for m_i in set(m) ])

try:

    # List of Small odd primes, L := [l_0, ..., l_{n-1}]
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

print("// All the experiments are assuming S = %1.3f x M and a = %1.3f x M.\n// The measures are given in millions of field operations.\n" % (SQR, ADD))
e = [1] * n
#MAX_iter = 16

sigma, kappa = 6, 1
SAMPLE = [ [0.0, 0.0, 0.0] ] * MAX_iter
bar = Bar('Experiments corresponding to SIMBA-' + str(sigma) + '-' + str(kappa) , max=MAX_iter)
B = list(A)
for main_i in range(MAX_iter):

    e = [ random.randint(0, m_i) for m_i in m]
    
    SET_FIELD_OPERATIONS_TO_ZERO()
    B = SIMBA(B, e, L, n, sigma, kappa, m, False)
    SAMPLE[main_i] = list(FIELD_OPERATIONS_PERFORMED)
    #print("Random(EllipticCurve(x^3 + 0x%X * x^2 + x)) * (p+1);" % coeff(b))
    
    bar.next()
bar.finish()

AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE ]), statistics.mean([ ops[1] for ops in SAMPLE ]), statistics.mean([ ops[2] for ops in SAMPLE ]) ]
print("Average number of field operations (GAE):\t%1.3f x M + %1.3f x S + %1.3f x a := %1.3f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )

tmp = list(set(m))
SAMPLE = [ [0.0, 0.0, 0.0] ] * MAX_iter
SAMPLE_FULL_TORSION_POINT = [ [0.0, 0.0, 0.0] ] * MAX_iter
bar = Bar('Experiments corresponding to this work' , max=MAX_iter)

for main_i in range(MAX_iter):

    e = [ random.randint(0, m_i) for m_i in m]
    
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
    
    #print("Random(EllipticCurve(x^3 + 0x%X * x^2 + x)) * (p+1);" % coeff(B))
    
    bar.next()
bar.finish()

# Fully constant-time group action evaluation (only one strategy evaluation is required)
if (len(tmp) == 1) and (tmp[0] == 1):

    SET_FIELD_OPERATIONS_TO_ZERO()
    B, m_tmp, e_tmp = evaluate_strategy(B, T_p, L_out[0], S_out[0], n, m, e)
    RUNNING_TIME = list(FIELD_OPERATIONS_PERFORMED)
    print("Fixed number of field operations (GAE):\t\t%1.3f x M + %1.3f x S + %1.3f x a := %1.3f x M" % (RUNNING_TIME[0] / (10.0**6), RUNNING_TIME[1] / (10.0**6), RUNNING_TIME[2] / (10.0**6), measure(RUNNING_TIME) / (10.0**6)) )
    
    AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE_FULL_TORSION_POINT ]), statistics.mean([ ops[1] for ops in SAMPLE_FULL_TORSION_POINT ]), statistics.mean([ ops[2] for ops in SAMPLE_FULL_TORSION_POINT ]) ]
    print("Randomly selection of full torsion points:\t%1.3f x M + %1.3f x S + %1.3f x a := %1.3f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )
else:
    AVERAGE = [statistics.mean([ ops[0] for ops in SAMPLE ]), statistics.mean([ ops[1] for ops in SAMPLE ]), statistics.mean([ ops[2] for ops in SAMPLE ]) ]
    print("Average number of field operations (GAE):\t%1.3f x M + %1.3f x S + %1.3f x a := %1.3f x M\n" % (AVERAGE[0] / (10.0**6), AVERAGE[1] / (10.0**6), AVERAGE[2] / (10.0**6), measure(AVERAGE) / (10.0**6)) )

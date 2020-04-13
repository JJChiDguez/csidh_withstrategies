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
    Number of degree-(l_i) isogeny constructions to be performed: m_i
    ------------------------------------------------------------------------------------- '''

# ==========================================================================

# __________________________________________________________________________
# ====== [MCR style, CSIDH-512] each m_i corresponds with the given in MCR18
#m = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5]

# ===== Suitable bounds for this work
#m = [15, 18, 20, 21, 21, 22, 22, 22, 22, 22, 22, 19, 20, 22, 23, 23, 23, 23, 23, 23, 23, 21, 23, 20, 16, 16, 16, 15, 14, 12, 13, 12, 11, 11, 10, 10, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3]

# __________________________________________________________________________
# ====== [OAYT style, CSIDH-512] each m_i corresponds with the given in OAYT19
#m = [5, 6, 7, 7, 7, 7, 7, 8, 8, 8, 9, 10, 10, 10, 10, 9, 9, 9, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1]

# ===== Suitable bounds for this work
#m = [7, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 9, 11, 9, 8, 8, 8, 7, 7, 7, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]

# __________________________________________________________________________
# ====== [dummy-free style, CSIDH-512] each m_i corresponds with the given in MCR18 (it is the same as MCR style)
#m = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5]

# ===== Suitable bounds for this work
#m = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5]

m = [15, 18, 20, 21, 21, 22, 22, 22, 22, 22, 22, 19, 20, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 19, 16, 16, 16, 15, 14, 12, 13, 12, 11, 11, 11, 9, 9, 9, 9, 8, 8, 8, 8, 7, 8, 6, 6, 6, 6, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3]

m.reverse()
# ---
k = 3   # Number of rows (format of the list)
# ---

C_m, L_m, R_m, S_m, t_m = strategy_block_cost(L[::-1], m)

L_MCR = [349, 347, 337, 331, 317, 313, 311, 307, 
293, 283, 281, 277, 271, 269, 263, 257, 
251, 241, 239, 233, 229, 227, 223, 211, 
199, 197, 193, 191, 181, 179, 173, 167, 
163, 157, 151, 149, 139, 137, 131, 127, 
113, 109, 107, 103, 101,  97,  89,  83,
 79,  73,  71,  67,  61,  59,  53,  47,
 43,  41,  37,  31,  29,  23,  19,  17,
 13,  11,   7,   5,   3, 587, 373, 367, 
359, 353]

print("#ifndef _STRATEGIES_H_")
print("#define _STRATEGIES_H_")
print("")
for i in range(len(t_m)):
    print("// -----------------------------------------------------------------------------------------------------------------------------------");
    print("// Strategy number %d\n" % (i))
    L_string = "static uint8_t L%d[] " % (i)
    R_string = "static uint8_t W%d[] " % (i)
    S_string = "static uint8_t S%d[] " % (i)
    
    printl(L_string, [ L_MCR.index(l_i) for l_i in L_m[i]][::-1], len(L_m[i]) // k + 1)
    if(R_m[i] != []):
        printl(R_string, [ L_MCR.index(r_i) for r_i in R_m[i]][::-1], len(R_m[i]) // k + 1)
    else:
        print(R_string +  " = {};")
    if(S_m[i] != []):
        printl(S_string, S_m[i], len(S_m[i]) // k + 1)
    else:
        print(S_string +  " = {};")
    

print("\n")
print("// -----------------------------------------------------------------------------------------------------------------------------------");
print("// -----------------------------------------------------------------------------------------------------------------------------------");
print("#define NUMBER_OF_DIFFERENT_STRATEGIES  %d" % len(L_m))
print("")
L_string = "static uint8_t *L_STRATEGY[NUMBER_OF_DIFFERENT_STRATEGIES] = {\n\t"
R_string = "static uint8_t *W_STRATEGY[NUMBER_OF_DIFFERENT_STRATEGIES] = {\n\t"
S_string = "static uint8_t *S[NUMBER_OF_DIFFERENT_STRATEGIES] = {\n\t"
for i in range(len(L_m) - 1):
    L_string = L_string + "L%d, " % (i)
    R_string = R_string + "W%d, " % (i)
    S_string = S_string + "S%d, " % (i)

L_string = L_string + "L%d\n\t};" % (len(L_m) - 1)
R_string = R_string + "W%d\n\t};" % (len(L_m) - 1)
S_string = S_string + "S%d\n\t};" % (len(L_m) - 1)

print("// L_STRATEGY[i] determines the small odd primes l_i per each strategy")
print(L_string)
print("\n// W_STRATEGY[i] determines L \ L_STRATEGY[i]")
print(R_string)
print("\n// S_STRATEGY[i] determines the optimal strategy for L_STRATEGY[i]")
print(S_string)

print("\n// Number of primes for each strategy")
printl("static uint8_t NUMBER_OF_PRIMES[]", [ len(L_m[i]) for i in range(len(L_m)) ], len(L_m) )
print("")
print("// Number of rounds per each different strategy")
printl("static uint8_t ROUNDS[]", t_m, len(t_m))

print("")
print("// Maximum number of degree-(l_i) isogeny constructions")
printl("static uint8_t B[]", [ m[L[::-1].index(l)] for l in L_MCR ], n // k + 1)

print("\n\n// Expected running-time ~ %f, and security %f.\n" % (measure(C_m), security(m, n)))
print("\n#endif /* all the required in the strategies for the approach with dummy operations */")

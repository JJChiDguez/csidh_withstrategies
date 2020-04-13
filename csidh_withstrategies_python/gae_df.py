from fp import *
from montgomery import *
from isogeny_hybrid import *

# random_key(m) implements an uniform random sample from S(m_1) x S(m_2) x ... x S(m_n)
random_key = lambda m: [ 2 * ( random.randint(0, m_i) - (m_i // 2) ) - (m_i % 2) for m_i in m]

'''
    security()
    inputs : the list M of maximum number of degree-(l_i) isogeny constructions 
             to be performed, and the length of M
    output : bits of security that M brings
'''
security = lambda M, n: sum(list([ log(M[i] + 1,2) for i in range(n)]));

# In order to achieve efficiency, the optimal strategies and their cost are saved in two global dictionaries (hash tables)
S = { 1: {} }    # Initialization of each strategy
C = { 1: {} }    # Initialization of the costs: 0.

for i in range(n):

    S[1][tuple([L[i]])] = [];           # Strategy with a list with only one element (a small odd prime number l_i)
    C[1][tuple([L[i]])] = C_xISOG[i];   # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]

for i in range(2, n+1):

    C[i] = { }
    S[i] = { }

'''
    dynamic_programming_algorithm():
    inputs: the list of small odd primes to be processed and its length
    output: the optimal strategy and its cost of the input list of small odd primes
'''
def dynamic_programming_algorithm(L, n):
    global S, C
    # If the approach uses dummy operations, to set DUMMY = 2.0;
    # otherwise, to set DUMMY = 1.0 (dummy free approach);

    if( len(L) != n ):

        # If the list of prime numbers doesn't have size n, then we return [],-1
        print("error:\tthe list of prime numbers has different size from %d." % n);
        return [], -1;
    else:

        # Assuming #L = n, we proceed.
        get_neighboring_sets = lambda L, k: [ tuple(L[i:i+k]) for i in range(n-k+1)] # This function computes all the k-tuple: (l_1, l_2, ..., l_{k)),
                                                                                     # (l_2, l_3, ..., l_{k+1)), ..., (l_{n-k}, l_{n-k+1, ..., l_{n)).
        for i in range(2, n+1):
 
            for Tuple in get_neighboring_sets(L, i):

                if C[i].get(Tuple) is None:

                    alpha = [ (b, 
                        C[len(Tuple[:b])][Tuple[:b]] +       # Subtriangle on the right side with b leaves
                        C[len(Tuple[b:])][Tuple[b:]] +       # Subtriangle on the left side with (i - b) leaves
                        2.0*sum([ C_xMUL[global_L.index(t)] for t in Tuple[:b] ])   +   # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                        2.0*sum([ C_xEVAL[global_L.index(t)] for t in Tuple[b:] ])  +   # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                        1.0*sum([ C_xMUL[global_L.index(t)] for t in Tuple[b:] ])       # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                        ) for b in range(1, i - 1)
                        ] +\
                        [ (i - 1,
                        C[i - 1][Tuple[:(i - 1)]] +       # Subtriangle on the right side with (i - 1) leaves
                        C[1][Tuple[(i - 1):]]     +       # Subtriangle on the left side with 1 leaf (only one vertex)
                        1.0*sum([ C_xMUL[global_L.index(t)] for t in Tuple[:(i - 1)] ]) + # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with 1 leaf
                        2.0*C_xEVAL[global_L.index(Tuple[i - 1])]  +                      # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                        1.0*C_xMUL[global_L.index(Tuple[i - 1])]                          # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - 1) leaves
                        )
                        ]
                    b, C[i][Tuple] = min(alpha, key=lambda t: measure(t[1]))     # We save the minimal cost corresponding to the triangle with leaves Tuple
                    S[i][Tuple] = [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]    # We save the optimal strategy corresponding to the triangle with leaves Tuple

        return S[n][tuple(L)], C[n][tuple(L)]   # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to C_xISOG[global_L.index(L[0])].


'''
    evaluate_strategy():
    inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
             a list of two projective Montgomery x-coordinate torsion-(l_1 x ... l_n) points x(T_+) := XP/ZP
             and x(T_-) := XM/ZM, the list of small odd primes [l_1, l_2, ..., l_n], an strategy and length
             of the given list of small odd primes, maximum number of degree-l_i isogeny constructions, and
             the secret integer vector to be evaluated
    output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
             is E / <P>, the new maximum and current number of degree-l_i isogeny constructions to be performed
             after the strategy evaluation.
    
    NOTE: T_+ belongs to E[pi - 1] and T_- belongs to E[pi + 1]. In particular, P = [T_-, T_+]
'''
def evaluate_strategy(E, P, L, strategy, n, m, e):

    v = list(m)
    u = list(e)
    ramifications = []
    moves = [0]         # moves: this list determines whether an isogeny construction must be performed
    k = 0               # k: current element of the strategy
    t = 1               # sign to be flip
    
    ramifications.append(list(P))   # list of pair of points (T_-, T_+) such that T_- in E[\pi + 1] and T_+ in E[\pi - 1]
    E_i = list(E)
    for i in range(len(strategy)):

        pos = global_L.index(L[n - 1 - i])      # Current element of global_L to be required

        # Reaching the vertex (n - 1 - i, i)

        # Vertical edges
        prev = sum(moves)
        while (prev + strategy[k]) < (n - 1 - i):

            moves.append(strategy[k])                   # Number of vertical edges to be performed
            T = list( [list(ramifications[-1][0]), list(ramifications[-1][1])] )
            for j in range(prev, prev + strategy[k], 1):
                T = list([xMUL(T[0], E_i, global_L.index(L[j])), xMUL(T[1], E_i, global_L.index(L[j])) ])

            ramifications.append(list([list(T[0]), list(T[1])]))
            prev += strategy[k]
            k += 1

        # Vertical edges without ramifications
        s_i = sign(e[pos])                                  # Sign of e[pos]
        c_i = (s_i + 1) // 2                                # Constant-swap of T_+ and T_-
        
        if prev < (n - 1 - i):

            moves.append(strategy[k])                           # Number of vertical edges (without ramifications) to be performed            
            T = list( [list(ramifications[-1][0]), list(ramifications[-1][1])] )
            
            T[0][0], T[1][0] = fp_cswap(T[0][0], T[1][0], c_i)
            T[0][1], T[1][1] = fp_cswap(T[0][1], T[1][1], c_i)
            for j in range(prev, prev + strategy[k], 1):
                T[0] = list(xMUL(T[0], E_i, global_L.index(L[j])))    # A single scalar multiplication is required
            
            T[0][0], T[1][0] = fp_cswap(T[0][0], T[1][0], c_i)
            T[0][1], T[1][1] = fp_cswap(T[0][1], T[1][1], c_i)
            
            ramifications.append(list([list(T[0]), list(T[1])]))
            prev += strategy[k]
            k += 1

        # At this point, vertex (n - 1 - i, i) has been reached
        if v[pos] > 0:     # Maximum number of degree-l_{pos} isogeny constructions?

            # At this step, ramifications[-1] is the i-th leaf
            # Swap the torsion-(l_{n-j}) elliptic curve points on the current leaf
            #                           T_- <- ramifications[-1][0]
            #                           T_+ <- ramifications[-1][1]
            ramifications[-1][0][0], ramifications[-1][1][0] = fp_cswap(ramifications[-1][0][0], ramifications[-1][1][0], c_i)  # XT_+ <-> XT_-
            ramifications[-1][0][1], ramifications[-1][1][1] = fp_cswap(ramifications[-1][0][1], ramifications[-1][1][1], c_i)  # ZT_+ <-> ZT_-

            if isinfinity(ramifications[-1][0]) == False:

                K = KPs(ramifications[-1][0], E_i, pos)
                
                # New isogeny construction
                E_i = xISOG(E_i, K, pos)
                
                # Horizontal edges are performed
                for j in range(0, len(moves) - 1, 1):
                
                    # Swap points in E[\pi - 1] and E[\pi + 1]
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(ramifications[j][0][0], ramifications[j][1][0], c_i)
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(ramifications[j][0][1], ramifications[j][1][1], c_i)

                    # T_- (or T_)
                    ramifications[j][0] = xEVAL(K, ramifications[j][0], pos)
                    # T_+ or (T_+)
                    ramifications[j][1] = xEVAL(K, ramifications[j][1], pos)
                    
                    # Multiplying by the current small odd prime number in order to decrease the order of both points T_+ and T_-
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)

                    # Undo swap of points in E[\pi - 1] and E[\pi + 1]
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(ramifications[j][0][0], ramifications[j][1][0], c_i)
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(ramifications[j][0][1], ramifications[j][1][1], c_i)

                b_i = isequal[u[pos] == 1] ^ isequal[u[pos] == -1]  # 0 if u[pos] != +1,-1; otherwise, 1 [This ask must be performed in constant-time]
                v[pos] -= 1
                u[pos] -= (s_i * (1 + b_i))     # reduced by 1 unit if it was different from +1 and -1; otherwise, reduced by 2 units

            else:

                for j in range(0, len(moves) - 1, 1):

                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(ramifications[j][0][0], ramifications[j][1][0], c_i)
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(ramifications[j][0][1], ramifications[j][1][1], c_i)
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(ramifications[j][0][0], ramifications[j][1][0], c_i)
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(ramifications[j][0][1], ramifications[j][1][1], c_i)
        else:
            # This branch only depends on randomness

            # At this step, ramifications[-1] is the i-th leaf
            ramifications[-1][0][0], ramifications[-1][1][0] = fp_cswap(ramifications[-1][0][0], ramifications[-1][1][0], c_i)
            ramifications[-1][0][1], ramifications[-1][1][1] = fp_cswap(ramifications[-1][0][1], ramifications[-1][1][1], c_i)

            if isinfinity(ramifications[-1][0]) == False:

                for j in range(0, len(moves) - 1, 1):

                    ramifications[j][0] = xMUL(ramifications[j][0], E_i, pos)
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)
            else:

                for j in range(0, len(moves) - 1, 1):

                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(ramifications[j][0][0], ramifications[j][1][0], c_i)
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(ramifications[j][0][1], ramifications[j][1][1], c_i)
                    ramifications[j][1] = xMUL(ramifications[j][1], E_i, pos)
                    ramifications[j][0][0], ramifications[j][1][0] = fp_cswap(ramifications[j][0][0], ramifications[j][1][0], c_i)
                    ramifications[j][0][1], ramifications[j][1][1] = fp_cswap(ramifications[j][0][1], ramifications[j][1][1], c_i)

        moves.pop()
        ramifications.pop()

    pos = global_L.index(L[0])                                      # Current element of global_L to be required
    s_i = sign(e[pos])                            # Sign of e[pos]
    c_i = (s_i + 1) // 2                          # Constant-swap of T_+ and T_-

    # T_+ or T_- ?
    # Root
    ramifications[0][0][0], ramifications[0][1][0] = fp_cswap(ramifications[0][0][0], ramifications[0][1][0], c_i)
    ramifications[0][0][1], ramifications[0][1][1] = fp_cswap(ramifications[0][0][1], ramifications[0][1][1], c_i)
    
    if isinfinity(ramifications[0][0]) == False:

        if m[pos] > 0:
        
            K = KPs(ramifications[0][0], E_i, pos)
            E_i = xISOG(E_i, K, pos)

            b_i = isequal[u[pos] == 1] ^ isequal[u[pos] == -1]  # 0 if u[pos] != +1,-1; otherwise, 1 [This ask must be performed in constant-time]
            v[pos] -= 1
            u[pos] -= (s_i * (1 + b_i))     # reduced by 1 unit if it was different from +1 and -1; otherwise, reduced by 2 units
        
    return E_i, v, u

# Next function computes the batches required in SIMBA
BATCHES = lambda L, n, sigma: [ [ L[j] for j in range(i, n, sigma) ] for i in range(sigma) ]

'''
    SIMBA():
    inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
             the secret integer vector to be evaluated, the list of small odd primes [l_1, l_2, ..., l_n] 
             and its length, number of batches and their number of times to be used, the maximum number of 
             degree-l_i isogeny constructions, and a bool variable that determines if optimal strategies are
             going to be used
    output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
             corresponds with the image of the group action evaluation.
'''
def SIMBA(A, e, L, n, sigma, kappa, m, OPTIMAL):

    E_k = list(A)

    # Batches to be used
    batches = BATCHES(L, n, sigma)
    # Complement of each batch
    remainders = [ [l for l in L if l not in batch] for batch in batches ]

    for j in range(kappa):

        for i in range(sigma):

            batch_size = len(batches[i])
            if batch_size > 0:
                T_p, T_m = elligator(E_k)
                T_p = xDBL(T_p, E_k)
                T_p = xDBL(T_p, E_k)
                T_m = xDBL(T_m, E_k)
                T_m = xDBL(T_m, E_k)

                for l in remainders[i]:
                    T_p = xMUL(T_p, E_k, global_L.index(l))
                    T_m = xMUL(T_m, E_k, global_L.index(l))

                if OPTIMAL == True:
                    St, Ct = dynamic_programming_algorithm(batches[i], batch_size)  # Optimal strategy
                else:
                    St = list(range(batch_size - 1, 0, -1))                         # Multiplicative strategy

                E_k, m, e = evaluate_strategy(E_k, list([list(T_m), list(T_p)]), batches[i], St, batch_size, m, e)

                # If the maximum of degree-(l_k) has been reached then the current batch (and its complement) must be updated
                tmp_batch     = [ batches[i][k] for k in range(batch_size) if m[global_L.index(batches[i][k])] > 0 ]
                tmp_remainder = [ batches[i][k] for k in range(batch_size) if m[global_L.index(batches[i][k])] == 0]

                batches[i] = list(tmp_batch)                    # Removing elements from the batch
                remainders[i] = list(remainders[i] + tmp_remainder)   # Adding elements to the complement of the batch

    # Multiplicative strategy on the set of unreached small odd prime numbers
    unreached_sop = [ global_L[i] for i in range(len(global_L)) if m[i] > 0 ]
    remainder_sop = [ l for l in global_L if l not in unreached_sop ]
    #print(unreached_sop)
    #print(remainder_sop)
    while len(unreached_sop) > 0:
        #print(e)

        T_p, T_m = elligator(E_k)
        T_p = xDBL(T_p, E_k)
        T_p = xDBL(T_p, E_k)
        T_m = xDBL(T_m, E_k)
        T_m = xDBL(T_m, E_k)

        for l in remainder_sop:
                T_p = xMUL(T_p, E_k, global_L.index(l))
                T_m = xMUL(T_m, E_k, global_L.index(l))

        current_n = len(unreached_sop)
        E_k, m, e = evaluate_strategy(E_k, list([list(T_m), list(T_p)]), unreached_sop, list(range(current_n - 1, 0, -1)), current_n, m, e)

        # If the maximum of degree-(l_k) has been reached then the current batch (and its complement) must be updated
        tmp_unreached = [ unreached_sop[k] for k in range(current_n) if m[global_L.index(unreached_sop[k])] > 0 ]
        tmp_remainder = [ unreached_sop[k] for k in range(current_n) if m[global_L.index(unreached_sop[k])] == 0]

        unreached_sop = list( tmp_unreached )                    # Removing elements from the batch
        remainder_sop = remainder_sop + tmp_remainder   # Adding elements to the complement of the batch

    return E_k

'''
    geometric_serie()
    inputs: and integer m, and a prime number l
    output: the nearest integer to
                  l
            m x -----
                l - 1
'''
def geometric_serie(m, l):
    
    l_float = float(l)
    m_float = float(m);
    return floor( (m_float * l_float) / (l_float - 1.0) + 0.5 )

'''
    filtered()
    inputs : a list L and a sublist SL of L
    output : L \ SL
'''
filtered = lambda List, sublist: [e for e in List if e not in sublist]

'''
    rounds()
    inputs : an integer vector (maximum number of isogeny constructions to be performed), 
             and the length of the vector
    output : the subset of (indexes of) the small odd primes that determines the optimal 
             strategy to be used, the number of times that each strategy will be used, and
             the complement of each subset (with respect to the set of all the small odd primes)
'''
def rounds(e, n):
    tmp_N = range(n)
    tmp_e = list(e)
    rounds_out = []
    sublists_L = []
    sublists_C = []
    while [e_i for e_i in tmp_e if e_i > 0] != []:
        e_min = min([e_i for e_i in tmp_e if e_i > 0])
        rounds_out.append(e_min)
        sublists_L.append([i for i in tmp_N if tmp_e[i] >= e_min])
        sublists_C.append(filtered(tmp_N, sublists_L[len(sublists_L) - 1]))
        tmp_e = [ (tmp_e[i] - e_min) for i in tmp_N ]
    
    return rounds_out, sublists_L, sublists_C;

'''
    GAE():
    inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
             the secret integer vector to be evaluated, all the sublist of small odd primes to be required and 
             their complement of each sublist, the (optimal) strategies corresponding to each sublist and their 
             number of times to be evaluated, and the maximum number of degree-l_i isogeny constructions
    output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
             corresponds with the image of the group action evaluation.

    NOTE: GAE comes from Group Action Evaluation, and the input sublists are determined by the output of function
          rounds().
          THIS IS THE IMPLEMENTATION OF OUR PROPOSED STRATEGY METHOD.
'''
def GAE(A, e, L, R, St, r, m):

    E_k = list(A)
    n = len(L)

    for j in range(0, n, 1):

        for k in range(0, r[j], 1):

            T_p, T_m = elligator(E_k)
            T_p = xDBL(T_p, E_k)
            T_p = xDBL(T_p, E_k)
            T_m = xDBL(T_m, E_k)
            T_m = xDBL(T_m, E_k)

            for l in R[j]:
                T_p = xMUL(T_p, E_k, global_L.index(l))
                T_m = xMUL(T_m, E_k, global_L.index(l))

            E_k, m, e = evaluate_strategy(E_k, list([list(T_m), list(T_p)]), L[j], St[j], len(L[j]), m, e)

    # Multiplicative strategy on the set of unreached small odd prime numbers
    unreached_sop = [ global_L[i] for i in range(len(global_L)) if m[i] > 0 ]
    remainder_sop = [ l for l in global_L if l not in unreached_sop ]

    while len(unreached_sop) > 0:

        T_p, T_m = elligator(E_k)
        T_p = xDBL(T_p, E_k)
        T_p = xDBL(T_p, E_k)
        T_m = xDBL(T_m, E_k)
        T_m = xDBL(T_m, E_k)

        for l in remainder_sop:
                T_p = xMUL(T_p, E_k, global_L.index(l))
                T_m = xMUL(T_m, E_k, global_L.index(l))

        current_n = len(unreached_sop)
        E_k, m, e = evaluate_strategy(E_k, list([list(T_m), list(T_p)]), unreached_sop, list(range(current_n - 1, 0, -1)), current_n, m, e)

        # If the maximum of degree-(l_k) has been reached then the current batch (and its complement) must be updated
        tmp_unreached = [ unreached_sop[k] for k in range(current_n) if m[global_L.index(unreached_sop[k])] > 0 ]
        tmp_remainder = [ unreached_sop[k] for k in range(current_n) if m[global_L.index(unreached_sop[k])] == 0]

        unreached_sop = list( tmp_unreached )                    # Removing elements from the batch
        remainder_sop = remainder_sop + tmp_remainder   # Adding elements to the complement of the batch

    return E_k


######################################################################################################################
# Next functions are used for computing optimal bounds
basis = np.eye(n, dtype = int)

# Next function computes the expected cost of our approach by assuming we have full torsion points
def strategy_block_cost(L, e):

    elligator_cost = np.array([7.0, 3.0, 10.0])    # Elligator cost
    mul_fp_by_four = np.array([8.0, 4.0, 8.0])     # Cost of computing x([4]P)

    n = len(L)
    e_prime = [ geometric_serie(e[k], L[k]) for k in range(n) ]

    tmp_r, tmp_Ls, tmp_Cs = rounds(e_prime, n)

    C_e = np.array([0.0, 0.0, 0.0])
    S_out = []
    L_out = []
    R_out = []
    for j in range(len(tmp_r)):

        R_out.append([L[k] for k in tmp_Cs[j]])
        L_out.append([L[k] for k in tmp_Ls[j]])

        bo_C = 2.0*sum([ C_xMUL[global_L.index(L[k])] for k in tmp_Cs[j] ])
        S_tmp, go_C = dynamic_programming_algorithm([ L[k] for k in tmp_Ls[j]] , len(tmp_Ls[j]))

        S_out.append(S_tmp)
        C_e += (go_C + bo_C + elligator_cost + 2.0*mul_fp_by_four) * tmp_r[j]

    return C_e, L_out, R_out, S_out, tmp_r;
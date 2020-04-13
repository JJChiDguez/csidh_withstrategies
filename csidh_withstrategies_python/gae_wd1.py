from fp import *
from montgomery import *
from isogeny_hybrid import *

random_key = lambda m: [ random.randint(0, m_i) for m_i in m]

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

    S[1][tuple([L[i]])] = [];                                                       # Strategy with a list with only one element (a small odd prime number l_i)
    C[1][tuple([L[i]])] = C_xISOG[i] - C_xMUL[i] + 2.0*np.array([4.0, 2.0, 6.0])    # For catching the weigth of horizontal edges of the form [(0,j),(0,j+1)]

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
                        C[len(Tuple[:b])][Tuple[:b]] +       # Subtriangle on the left side with b leaves
                        C[len(Tuple[b:])][Tuple[b:]] +       # Subtriangle on the right side with (i - b) leaves
                        1.0*sum([ C_xMUL[global_L.index(t)] for t in Tuple[:b] ])   +   # Weights corresponding with vertical edges required for connecting the vertex (0,0) with the subtriangle with b leaves
                        1.0*sum([ C_xEVAL[global_L.index(t)] for t in Tuple[b:] ])  +   # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                        1.0*sum([ C_xMUL[global_L.index(t)] for t in Tuple[b:] ])       # Weights corresponding with horizontal edges required for connecting the vertex (0,0) with the subtriangle with (i - b) leaves
                        ) for b in range(1, i)
                        ]
                    b, C[i][Tuple] = min(alpha, key=lambda t: measure(t[1]))     # We save the minimal cost corresponding to the triangle with leaves Tuple
                    S[i][Tuple] = [b] + S[i - b][Tuple[b:]] + S[b][Tuple[:b]]    # We save the optimal strategy corresponding to the triangle with leaves Tuple

        return S[n][tuple(L)], C[n][tuple(L)] + C_xMUL[global_L.index(L[0])] - 2.0*np.array([4.0, 2.0, 6.0])  # The weight of the horizontal edges [(0,n-1),(0,n)] must be equal to C_xISOG[global_L.index(L[0])].


'''
    evaluate_strategy():
    inputs : a projective Montgomery constants A24:= A + 2C and C24:=4C where E : y^2 = x^3 + (A/C)*x^2 + x,
             a projective Montgomery x-coordinate torsion-(l_1 x ... l_n) point x(P) := XP/ZP, the list of
             small odd primes [l_1, l_2, ..., l_n], an strategy and length of the given list of small odd 
             primes, maximum number of degree-l_i isogeny constructions, and the secret integer vector to be
             evaluated
    output : the projective Montgomery constants a24:= a + 2c and c24:=4c where E': y^2 = x^3 + (a/c)*x^2 + x
             is E / <P>, the new maximum and current number of degree-l_i isogeny constructions to be performed
             after the strategy evaluation.
'''
def evaluate_strategy(E, P, L, strategy, n, m, e):

    v = list(m)
    u = list(e)
    ramifications = []
    moves = [0]         # moves: this list determines whether an isogeny construction must be performed
    k = 0               # k: current element of the strategy

    ramifications.append(P)
    E_i = list(E)
    for i in range(len(strategy)):

        pos = global_L.index(L[n - 1 - i])      # Current element of global_L to be required

        # Reaching the vertex (n - 1 - i, i)

        # Vertical edges
        prev = sum(moves)
        while prev < (n - 1 - i):

            moves.append(strategy[k])                   # Number of vertical edges to be performed
            T = list( ramifications[-1] )   # New ramification
            for j in range(prev, prev + strategy[k], 1):
                T = xMUL(T, E_i, global_L.index(L[j]))

            ramifications.append(T)
            prev += strategy[k]
            k += 1

        # At this point, vertex (n - 1 - i, i) has been reached
        if v[pos] > 0:     # Maximum number of degree-l_{pos} isogeny constructions?

            # At this step, ramifications[-1] is the i-th leaf
            if isinfinity(ramifications[-1]) == False:

                # Dummy or NOT Dummy degree-(l_{n-1-i}) isogeny construction, that's the question?
                b_i = isequal[u[pos] == 0]

                ramifications[-1][0], ramifications[0][0] = fp_cswap(ramifications[-1][0], ramifications[0][0], b_i)
                ramifications[-1][1], ramifications[0][1] = fp_cswap(ramifications[-1][1], ramifications[0][1], b_i)

                K = KPs(ramifications[-1], E_i, pos)

                ramifications[-1][0], ramifications[0][0] = fp_cswap(ramifications[-1][0], ramifications[0][0], b_i)
                ramifications[-1][1], ramifications[0][1] = fp_cswap(ramifications[-1][1], ramifications[0][1], b_i)

                # New isogeny construction
                C_i = xISOG(E_i, K, pos)

                # Next, the horizontal edge [(0,i),(0,i+1)] is performed
                d_i = (global_L[pos] - 1) // 2
                mask = isequal[(global_L[pos] == 3)]                          # catching special case when l = 3

                Z = yADD(K[(d_i + mask) - 1], K[0], K[(d_i + mask) - 2])      # y([d_i + 1]K[0])
                Z[0], K[d_i][0] = fp_cswap(Z[0], K[d_i][0], mask ^ 1)
                Z[1], K[d_i][1] = fp_cswap(Z[1], K[d_i][1], mask ^ 1)

                T = yADD(K[d_i], K[d_i - 1], K[0])                  # y([2*d_i + 1]K[0]) := y([l_i]K[0])
                T = [ fp_add(T[1], T[0]), fp_sub(T[1], T[0]) ]      # x([l_i]K[0])

                ramifications[0] = xEVAL(K, ramifications[0], pos)
                T[0], ramifications[0][0] = fp_cswap(T[0], ramifications[0][0], b_i)
                T[1], ramifications[0][1] = fp_cswap(T[1], ramifications[0][1], b_i)

                # The remainder horizontal edges are performed
                for j in range(1, len(moves) - 1, 1):

                    T = xMUL(ramifications[j], E_i, pos)
                    ramifications[j] = xEVAL(K, ramifications[j], pos)
                    T[0], ramifications[j][0] = fp_cswap(T[0], ramifications[j][0], b_i)
                    T[1], ramifications[j][1] = fp_cswap(T[1], ramifications[j][1], b_i)

                C_i[0], E_i[0] = fp_cswap(C_i[0], E_i[0], b_i ^ 1)
                C_i[1], E_i[1] = fp_cswap(C_i[1], E_i[1], b_i ^ 1)

                v[pos] -= 1
                u[pos] -= (b_i ^ 1)
        else:
            for j in range(0, len(moves) - 1, 1):

                ramifications[j] = xMUL(ramifications[j], E_i, pos)

        moves.pop()
        ramifications.pop()

    pos = global_L.index(L[0])                                      # Current element of global_L to be required
    if isinfinity(ramifications[0]) == False:

        if m[pos] > 0:

            b_i = isequal[e[pos] == 0]

            K = KPs(ramifications[0], E_i, pos)
            C_i = xISOG(E_i, K, pos)
            C_i[0], E_i[0] = fp_cswap(C_i[0], E_i[0], b_i ^ 1)
            C_i[1], E_i[1] = fp_cswap(C_i[1], E_i[1], b_i ^ 1)

            v[pos] -= 1
            u[pos] -= (b_i ^ 1)

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

                for l in remainders[i]:
                    T_p = xMUL(T_p, E_k, global_L.index(l))

                if OPTIMAL == True:
                    St, Ct = dynamic_programming_algorithm(batches[i], batch_size)  # Optimal strategy
                else:
                    St = list(range(batch_size - 1, 0, -1))                         # Multiplicative strategy

                E_k, m, e = evaluate_strategy(E_k, T_p, batches[i], St, batch_size, m, e)

                # If the maximum of degree-(l_k) has been reached then the current batch (and its complement) must be updated
                tmp_batch     = [ batches[i][k] for k in range(batch_size) if m[global_L.index(batches[i][k])] > 0 ]
                tmp_remainder = [ batches[i][k] for k in range(batch_size) if m[global_L.index(batches[i][k])] == 0]

                batches[i] = list(tmp_batch)                    # Removing elements from the batch
                remainders[i] = list(remainders[i] + tmp_remainder)   # Adding elements to the complement of the batch

    # Multiplicative strategy on the set of unreached small odd prime numbers
    unreached_sop = [ global_L[i] for i in range(len(global_L)) if m[i] > 0 ]
    remainder_sop = [ l for l in global_L if l not in unreached_sop ]
    while len(unreached_sop) > 0:

        T_p, T_m = elligator(E_k)
        T_p = xDBL(T_p, E_k)
        T_p = xDBL(T_p, E_k)

        for l in remainder_sop:
                T_p = xMUL(T_p, E_k, global_L.index(l))

        current_n = len(unreached_sop)
        E_k, m, e = evaluate_strategy(E_k, T_p, unreached_sop, list(range(current_n - 1, 0, -1)), current_n, m, e)

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

            for l in R[j]:
                T_p = xMUL(T_p, E_k, global_L.index(l))

            E_k, m, e = evaluate_strategy(E_k, T_p, L[j], St[j], len(L[j]), m, e)

    # Multiplicative strategy on the set of unreached small odd prime numbers
    unreached_sop = [ global_L[i] for i in range(len(global_L)) if m[i] > 0 ]
    remainder_sop = [ l for l in global_L if l not in unreached_sop ]

    while len(unreached_sop) > 0:

        T_p, T_m = elligator(E_k)
        T_p = xDBL(T_p, E_k)
        T_p = xDBL(T_p, E_k)

        for l in remainder_sop:
                T_p = xMUL(T_p, E_k, global_L.index(l))

        current_n = len(unreached_sop)
        E_k, m, e = evaluate_strategy(E_k, T_p, unreached_sop, list(range(current_n - 1, 0, -1)), current_n, m, e)

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

        bo_C = 1.0*sum([ C_xMUL[global_L.index(L[k])] for k in tmp_Cs[j] ])
        S_tmp, go_C = dynamic_programming_algorithm([ L[k] for k in tmp_Ls[j]] , len(tmp_Ls[j]))

        S_out.append(S_tmp)
        C_e += (go_C + bo_C + elligator_cost + 1.0*mul_fp_by_four) * tmp_r[j]

    return C_e, L_out, R_out, S_out, tmp_r;


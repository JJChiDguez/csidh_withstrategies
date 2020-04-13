from fp import *
from montgomery import *

''' ----------------------------------------------------------------------------
    KPs()
    input : the projective Montgomery x-coordinate points x(P) := XP/ZP, the 
            projective Montgomery constants A24:= A + 2C and C24:=4C where 
            E : y^2 = x^3 + (A/C)*x^2 + x, and a positive integer 0 <= i < n
    output: the list of projective Montgomery x-coordinate points x(P), x([2]P),
            x([3]P), ..., and x([d_i]P) where l_i = 2 * d_i + 1
    ---------------------------------------------------------------------------- '''
def KPs(P, A, i):

    d = (L[i] - 1) // 2

    K = [ [0, 0] for j in range(d + 1) ]
    K[0] = list(P)
    K[1] = xDBL(P, A)                           # 4M + 2S + 4a

    for j in range(2, d, 1):
        K[j] = xADD(K[j-1], K[0], K[j - 2])     # 4M + 2S + 6a

    return K                                    # 2(l - 3)M + (l - 3)S + 3(l - 3)a - 2a

''' --------------------------------------------------------------------------
    xISOG()
    input : the projective Montgomery constants A24:= A + 2C and C24:=4C where
            E : y^2 = x^3 + (A/C)*x^2 + x, the list of projective Montgomery
            x-coordinate points x(P), x([2]P), x([3]P), ..., and x([d_i]P)
            where l_i = 2 * d_i + 1, and a positive integer 0 <= i < n
    output: the projective Montgomery constants a24:= a + 2c and c24:=4c where
            E': y^2 = x^3 + (a/c)*x^2 + x is a degree-(l_i) isogenous curve to
            E, and the list of projective Twisted Edwards y-coordinate points
            y(P), y([2]P), y([3]P), ..., and y([d_i]P) where l_i = 2 * d_i + 1
    -------------------------------------------------------------------------- '''
def xISOG(A, K, i):

    d = (L[i] - 1) // 2

    pi = [1, 1]
    sg = [0, 1]

    S = [ [0,0] for j in range(d) ]

    for j in range(d):

        t = fp_mul(K[j][0], K[j][1])
        s_1 = fp_add(K[j][0], K[j][1]);
        s_2 = fp_sub(K[j][0], K[j][1]);
        S[j] = [s_1, s_2]

        sg_0 = fp_mul(s_1, s_2)
        sg_0 = fp_mul(sg_0, sg[1])
        tmp  = fp_mul(sg[0], t)

        sg = [ fp_add(sg_0, tmp), fp_mul(sg[1], t)]
        pi = [ fp_mul(pi[0], K[j][0]), fp_mul(pi[1], K[j][1])]

    A24 = fp_add(A[0], A[0])
    A24 = fp_sub(A24, A[1])
    A24 = fp_add(A24, A24)

    pi_0_squared = fp_sqr(pi[0])
    tmp_1 = fp_mul(A24, sg[1])
    tmp_2 = fp_mul(sg[0], A[1])
    tmp_3 = fp_add(tmp_2, tmp_2)
    tmp_2 = fp_add(tmp_2, tmp_3)
    tmp_2 = fp_add(tmp_2, tmp_2)
    tmp_1 = fp_sub(tmp_1, tmp_2)
    A24 = fp_mul(pi_0_squared, tmp_1)

    pi_1_squared = fp_sqr(pi[1])
    C24 = fp_mul(A[1], sg[1])
    C24 = fp_mul(pi_1_squared, C24)

    C24 = fp_add(C24, C24)
    A24 = fp_add(A24, C24)
    C24 = fp_add(C24, C24)

    return [A24, C24], S;   #end function;   // (7d + 5)M + 2S + (3d + 10)a

'''
xEVAL := function(P, K, d)
    Q := [P[1] + P[2], P[1] - P[2]];
    R := CrissCross(K[1][1], K[1][2], Q[1], Q[2]);
    for j:=2 to d do
        T := CrissCross(K[j][1], K[j][2], Q[1], Q[2]);
        R := [T[1] * R[1], T[2] * R[2]];
    end for;
    return [P[1] * (R[1]^2), P[2] * (R[2]^2)];
end function;   // (4d)M + 2S + 2(d + 1)a

'''
def xEVAL(K, P, i):

    d = (L[i] - 1) // 2            # Here, l = 2d + 1

    Q0 = fp_add(P[0], P[1])
    Q1 = fp_sub(P[0], P[1])
    R0, R1 = CrissCross(K[0][0], K[0][1], Q0, Q1)
    for j in range(1, d, 1):

        T0, T1 = CrissCross(K[j][0], K[j][1], Q0, Q1)
        R0 = fp_mul(T0, R0)
        R1 = fp_mul(T1, R1)

    R0 = fp_sqr(R0)
    R1 = fp_sqr(R1)
    X = fp_mul(P[0], R0)
    Z = fp_mul(P[1], R1)
    return [X, Z]                               # 2(l - 1)M + 2S + (l + 1)a

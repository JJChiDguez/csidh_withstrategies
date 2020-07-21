# Optimal strategies for CSIDH

This is joined work with:

        * Jesús-Javier Chi-Domínguez, and
        * Francisco Rodríguez-Henríquez
    
        This work uses the fact that the action used in the CSIDH protocol 
        can be computed with
                        1. shortest differential addition chains,
                        2. only Edwards curves (no hybrid between Montgomery and Edwards),
                        3. A projective elligator with u randomly chosen from {2, (p-1)/2}, and
                        4. the only two ``if statements'' used are for asking if a ``public'' 
                        point is the infinity point.
        The implemented methods are Simba with one and two torsion points,
        which require dummy operations.
        In addition, a dummy-free CSIDH approach was implemented where each 
        secret exponent e_i corresponding to l_i is randomly chosen from
                        1. the set of all odd integer number in {-b_i, ..., b_i} if b_i is odd, or
                        2. the set of all even integer number in {-b_i, ..., b_i} if b_i is even.
        Here, b_i is the number of degree-(l_i) isogenies to be computed.
    
        [Notes]
        [1]
                The Simba method implemented is based on the work of Michael Meyer, 
                Fabio Campos, Steffen Reith: 
                "On Lions and Elligators: An efficient constant-time implementation of CSIDH". 
                IACR Cryptology ePrint Archive 2018: 1198 (2018).
        [2]
                This implementation re-use the field arithmetic implemented in the
                original CSIDH paper by Wouter Castryck, Tanja Lange, Chloe Martindale, 
                Lorenz Panny, Joost Renes: 
                "CSIDH: An Efficient Post-Quantum Commutative Group Action". 
                ASIACRYPT (3) 2018: 395-427.
                (its eprint version can be download at https://eprint.iacr.org/2018/383)
        [3]
                This implementation re-use the given one in the following paper:
                Daniel Cervantes-Vázquez, Mathilde Chenu, Jesús-Javier. Chi-Domınguez, 
                Luca De Feo, Francisco Rodrıguez-Henrıquez, and Benjamin Smith, 
                “Stronger and Faster Side-Channel Protections for CSIDH”, 
                Progress in Cryptology - LATINCRYPT 2019. LNCS 11774 (2019), 173-193



This C code implementation was performed by:

        * Jesús-Javier Chi-Domínguez <jesus.chidominguez@tuni.fi, chidoys@gmail.com, jjchi@computacion.cs.cinvestav.mx>, and
        * Francisco Rodríguez-Henríquez <francisco@cs.cinvestav.mx>.

# C code

To compile the files you can do the following. First, you can use any version of gcc compiler (just set the variable CC as an input of the Makefile [variable CC is optional, gcc is set by default]).

## Testing a CSIDH protocol (key exchange protocol)

[Compilation]

```bash
    [SIMBA method]
            (Using dummy operations and one torsion point)
                    make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=SIMBA
            (Using dummy operations and two torsion points)
                    make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=SIMBA
            (Dummy-free approach and using two torsion points)
                    make csidh BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=SIMBA

    [Our proposed strategy method]
            (Using dummy operations and one torsion point)
                    make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=STRATEGY
            (Using dummy operations and two torsion points)
                    make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=STRATEGY
            (Dummy-free approach and using two torsion points)
                    make csidh BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=STRATEGY
```

[Execution]

```bash
./bin/csidh
```

## Running-time: number of field operations

[Compilation]

```bash
    [SIMBA method]
            (Using dummy operations and one torsion point)
                    make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=SIMBA
            (Using dummy operations and two torsion points)
                    make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=SIMBA
            (Dummy-free approach and using two torsion points)
                    make action_cost BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=SIMBA

    [Our proposed strategy method]
            (Using dummy operations and one torsion point)
                    make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=STRATEGY
            (Using dummy operations and two torsion points)
                    make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=STRATEGY
            (Dummy-free approach and using two torsion points)
                    make action_cost BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=STRATEGY
```
[Execution]
                ./bin/action_cost

## Running-time: number of clock cycles

[Compilation]

```bash
    [SIMBA method]
            (Using dummy operations and one torsion point)
                    make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=SIMBA
            (Using dummy operations and two torsion points)
                    make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=SIMBA
            (Dummy-free approach and using two torsion points)
                    make action_timing BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=SIMBA

    [Our proposed strategy method]
            (Using dummy operations and one torsion point)
                    make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=STRATEGY
            (Using dummy operations and two torsion points)
                    make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=STRATEGY
            (Dummy-free approach and using two torsion points)
                    make action_timing BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=STRATEGY
```

[Execution]
                ./bin/action_timing

## Cleaning data

```bash
    make clean
```

## Funding

This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 804476).
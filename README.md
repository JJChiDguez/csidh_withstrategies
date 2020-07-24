# Optimal strategies for CSIDH

CSIDH protocol was original proposed in 

```text
This implementation re-use the field arithmetic implemented in the
original CSIDH paper by Wouter Castryck, Tanja Lange, Chloe Martindale, 
Lorenz Panny, Joost Renes: 
"CSIDH: An Efficient Post-Quantum Commutative Group Action". 
ASIACRYPT (3) 2018: 395-427.
(its eprint version can be download at https://eprint.iacr.org/2018/383)            
```

This project gives a C and Python3 -codes implementations of CSIDH protocol. This work uses the fact that the action used in the CSIDH protocol can be computed with

1. shortest differential addition chains,
2. only Edwards curves (no hybrid between Montgomery and Edwards),
3. s projective elligator with u randomly chosen from `{2, (p-1)/2}`, and
4. the only two `if statements` used are for asking if a `public`  point is the infinity point.

The C-code implementation is an extension re-use the given one in the following paper:

```text
Daniel Cervantes-Vázquez, Mathilde Chenu, Jesús-Javier. Chi-Domınguez, 
Luca De Feo, Francisco Rodríguez-Henríquez, and Benjamin Smith, 
“Stronger and Faster Side-Channel Protections for CSIDH”, 
Progress in Cryptology - LATINCRYPT 2019. LNCS 11774 (2019), 173-193
```

To be more precise, the provided implementation compares the `SIMBA` method  with our proposed use of strategy applied on

1. MCR-style (one torsion point and dummy isogeny constructions):

   ```text
   Michael Meyer, Fabio Campos, Steffen Reith,
   "On Lions and Elligators: An efficient constant-time implementation of CSIDH". 
   Post-Quantum Cryptography - PQCrypto 2019. LNCS 11505 (2019), 307--325
   ```

2. OAYT-style (two torsion points and dummy isogeny constructions); and

   ```text
   Hiroshi Onuki, Yusuke Aikawa, Tsutomu Yamazaki, and Tsuyoshi Takagi,
   "(Short Paper) A Faster Constant-Time Algorithm of CSIDH Keeping Two Points".
   Advances in Information and Computer Security - IWSEC 2019, LNCS 11689 (2019), 23--33.
   (its eprint version can be download at https://eprint.iacr.org/2019/353)
   ```

3. CCCDRS-style (two torsion points and dummy-free approach).

# Compiling

First, you can use any version of gcc compiler (just set the variable CC as an input of the Makefile [variable CC is optional, gcc is set by default]). To compile the c-code implementation you can do the following steps. 

### Testing a CSIDH protocol (key exchange protocol)

```bash
# SIMBA method (using dummy operations and one torsion point)
make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=SIMBA
# SIMBA method (using dummy operations and two torsion points)
make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=SIMBA
# SIMBA method (dummy-free approach and using two torsion points)
make csidh BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=SIMBA

# This work (using dummy operations and one torsion point)
make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=STRATEGY
# This work (using dummy operations and two torsion points)
make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=STRATEGY
# This work (dummy-free approach and using two torsion points)
make csidh BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=STRATEGY
```

Once the desired compilation has been performed, just run

```bash
./bin/csidh
```

### Running-time: number of field operations

```bash
# SIMBA method (using dummy operations and one torsion point)
make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=SIMBA
# SIMBA method (using dummy operations and two torsion points)
make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=SIMBA
# SIMBA method (dummy-free approach and using two torsion points)
make action_cost BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=SIMBA

# This work (using dummy operations and one torsion point)
make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=STRATEGY
# This work (using dummy operations and two torsion points)
make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=STRATEGY
# This work (dummy-free approach and using two torsion points)
make action_cost BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=STRATEGY
```

Once the desired compilation has been performed, just run

```bash
./bin/action_cost
```

### Running-time: number of clock cycles

```bash
# SIMBA method (using dummy operations and one torsion point)
make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=SIMBA
# SIMBA method (using dummy operations and two torsion points)
make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=SIMBA
# SIMBA method (dummy-free approach and using two torsion points)
make action_timing BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=SIMBA

# This work (using dummy operations and one torsion point)
make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1 APPROACH=STRATEGY
# This work (using dummy operations and two torsion points)
make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2 APPROACH=STRATEGY
# This work (dummy-free approach and using two torsion points)
make action_timing BITLENGTH_OF_P=512 TYPE=DUMMYFREE APPROACH=STRATEGY
```

Once the desired compilation has been performed, just run

```bash
./bin/action_timing
```

### Cleaning data

```bash
    make clean
```

### Remarks

The Python3-code implementation has its own [README.md](/csidh_withstrategies_python/README.md).

## Authors

1. **_Jesús-Javier Chi-Domínguez_** <jesus.chidominguez@tuni.fi>, <chidoys@gmail.com>, <jjchi@computacion.cs.cinvestav.mx>, and
2. **_Francisco Rodríguez-Henríquez_** <francisco@cs.cinvestav.mx>.

## License

This project is licensed under the GNU general public license - see the [LICENSE](LICENSE) file for details

## Funding

This project has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 804476).
## Hutchinson _et al._ octave-based code configuration

The files contained in this directory corresponds with given one by Hutchinson _et al._ and can be downloaded from their [repository](https://github.com/AaronHutchinson/CSIDH).

```bash
# Their Octave and Matlab based codes can be downloaded by running the following
git clone https://github.com/AaronHutchinson/CSIDH.git
```

Notice, the **Octave**-based code provided by Hutchinson _et al._ in https://eprint.iacr.org/2019/1121 requires the **_CVX_** package ([CVX: Matlab Software for Disciplined Convex Programming](http://cvxr.com/cvx/)), which is provided by **Matlab**. However, having a **Matlab** license is not alway possible  (which was our case); for that reason, we decided to install and configure the **_CVX_** package on **GNU Octave**. 

### Installing GNU Octave and the CVX package

The following steps (commands) correspond with an installation on a Ubuntu Linux distribution. Installing **GNU Octave** can be installed by running the following in a Linux console:

```bash
sudo apt-get install octave
sudo apt-get install liboctave-dev
```

Once **GNU Octave** has been installed, next step is to download the CVX package from the `cvxr` git repository. That repository requires  `libopenblas-dev` library (the directory `octave_cvx` can be created wherever you decide):

```bash
sudo apt-get install libopenblas-dev
mkdir octave_cvx && cd octave_cvx
git clone -b rework https://github.com/cvxr/CVX.git/
cd CVX
rmdir sedumi sdpt3
git clone https://github.com/sqlp/sdpt3.git
git clone https://github.com/sqlp/sedumi.git
```

At this level, we are ready for installing the  **_CVX_** package. Now, one needs to change the function `octave_config_info` by `__octave_config_info__` in the file `cvx_version.m` (line 50). In the next steps, the flag `-rebuild` is required if you had installed before any **_CVX_** package version.

```bash
octave cvx_compile.m -rebuild
cd sedumi/
octave install_sedumi.m -nopath -rebuild
cd ..
cd sdpt3/
octave sdpt3.m -nopath -rebuild
cd ..
octave cvx_setup.m -rebuild
```

At this step, if everything has been correctly installed then **_CVX_** package is now installed on your computer. Nevertheless, if one tries to run the **Octave**-based code of Hutchinson _et al._, an error will be showed (we require some extra steps for configuring the **_CVX_** package); this is because the files `./lib/@cvxprob/display.m`, `./lib/@cvxdual/display.m`, `./lib/@cvx/display.m`, `./lib/@cvxtuple/display.m`, and `./lib/@cvxcnst/display.m` are using the expression `long = ~isequal(get(0,'FormatSpacing'),'compact');` and `get(0,'FormatSpacing')` is not a valid syntax in **GNU Octave**. In order to fix it, one can change the expression `long = ~isequal(get(0,'FormatSpacing'),'compact');` by

```octave
[~, fmt] = format ();
long = ~isequal(fmt,'compact');
```

Notice the directory `octave_cvx` cannot be moved or deleted (the **_CVX_** package is there).  Now, one should be able in running the **Octave**-based code provided by Hutchinson _et al._

### Remarks

The above steps were successfully performed in Lubuntu 20.04 distro.
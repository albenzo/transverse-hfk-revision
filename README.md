# transverse-hfk-revision
This is a cleanup of the transverseHFK.c to compute the invariant defined in 
"Legendrian knots, transverse knots, and combinatorial Floer homology" by 
P. S. Ozsvath, Z. Szabo, and D. P. Thurston. The original code is from that 
linked in the paper "Transverse knots distinguished by Knot Floer Homology" 
by L. Ng, P. S. Ozsvath, and D. P. Thurston.

There is a python library included that allows calls to the original
null_homologous_D0Q and null_homologous_D1Q methods as well as the 
invariants directly in a Tk window.

## transverseHFK
### Usage
To compile transverseHFK program run
```
$ make
```
The program will calculate whether or not x^-, x^+, delta_1(x^-), and delta_1(x^+) 
are null-homologous for the supplied knot. 

The program usage is
```
$ transverseHFK -i <ArcIndex> -X <List of Xs> -O <List of Os>
```
where the three parameters specify a knot via its grid diagram. Currently the
two lists must be input in the form `[0,1,2,3,4,5,6,7,8,9]` with no spaces.

For a full list of options run `transverseHFK --help`

### Installation and removal
To install the program type 
```
$ make install
```
To remove type 
```
$ make uninstall
```
Which will copy and remove `transverseHFK` to `/usr/bin` respectively.
If you do not wish to (or cannot) install to this location you can 
install by copying the executable to somewhere on your `$PATH`.

### Sample Output
```
$ transverseHFK -i 10 -X [10,5,8,6,3,7,2,4,9,1] -O [7,9,3,4,5,1,6,10,2,8]
    LL is null-homologous
    UR is NOT null-homologous
    D1[LL] is null-homologous
    D1[UR] is NOT null-homologous
```

## transHFK
### Usage
To build and install the python library run
```
$ python setup.py install
```
or
```
$ make python-install
```

Then it can be used by
```
$ python
>>> from transHFK import *
```
which exposes the `transHFK` and `Tk_transHFK` classes as well as
the `null_homologous_D0Q` and `null_homologous_D1Q` methods.
Usage information can be found via `help(<class or method name>)`

### Gridlink
To integrate transHFK with [gridlink](https://www.math.uic.edu/~culler/gridlink) 
run the following command in the terminal

```
patch <path to gridlink/gridlink.py> gridlink.patch
```

Then after installing gridlink via setup.py as normal and installing the transHFK python
library there will be a menu that allows you to run the transverseHFK 
invariants.

Note: Gridlink uses the opposite convention as to whether horizontal strands are
below or above for crossings. When gridlink passes our program it silently converts
the grid to match our convention.

## Platform specific notes
Building on mac requires argp-standalone to be installed. In addition if while
compiling the python libraries on mac you encounter a `duplicate symbol` bug this
can usually be resolved by running 
```
$ CFLAGS="fcommon" python setup.py install
```

Building natively on Windows is currently not supported. It is recommended that if
you wish to use the program on Windows that you work either inside the 
Windows subsystem for linux if you are on Windows 10 or to use Cygwin otherwise.

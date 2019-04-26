# transverse-hfk-revision
This project provides a cleaned up version to compute the invariant defined in 
"Legendrian knots, transverse knots, and combinatorial Floer homology" by 
P. S. Ozsvath, Z. Szabo, and D. P. Thurston. The original code written by P. Ozsvath 
is based on the paper "Transverse knots distinguished by Knot Floer Homology" 
by L. Ng, P. S. Ozsvath, and D. P. Thurston. The original code can be found
[here](https://services.math.duke.edu/~ng/math/programs.html).

In addition we provide support for computing a family of invariants described
in the paper (place holder).

There is a python library included that allows calls to the original
null_homologous_D0Q and null_homologous_D1Q methods as well as the 
invariants directly in a Tk window or from the python interpreter.

## transverseHFK
### Usage
To compile transverseHFK program run
```
$ make
```
The program will calculate whether or not x^-, x^+, delta_1(x^-), and delta_1(x^+) 
are null-homologous for the supplied knot. 

When used as
```
$ transverseHFK -i <ArcIndex> -X <List of Xs> -O <List of Os>
```
and for the new invariant is
```
$ transverseHFK -i <arcIndex> -X <List of Xs> -O <List of Os> -n <# of sheets>
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

## tHFK
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
>>> from tHFK import *
```
which exposes the `tHFK` and `Tk_tHFK` classes as well as
the `null_homologous_D0Q` and `null_homologous_D1Q` methods.
Usage information can be found via `help(<class or method name>)`

### Gridlink
To integrate tHFK with gridlink run the following command in the terminal

```
patch <path to gridlink/gridlink.py> gridlink.patch
```

Then after installing gridlink via setup.py as normal and installing the tHFK python
library there will be a menu that allows you to run the transverseHFK 
invariants.

Note: Gridlink uses the opposite convention as to whether horizontal strands are
below or above for crossings. When gridlink passes our program it silently converts
the grid to match our convention.

## Platform specific notes
Building the command-line executable on mac requires argp-standalone to be installed.

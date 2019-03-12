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
make
```
The program will calculate whether or not x^-, x^+, delta_1(x^-), and delta_1(x^+) 
are null-homologous for the supplied knot. 

The program usage is
```
transverseHFK -i <ArcIndex> -X <List of Xs> -O <List of Os>
```
where the three parameters specify a knot via its grid diagram. Currently the
two lists must be input in the form `[0,1,2,3,4,5,6,7,8,9]` with no spaces.

For a full list of options run `transverseHFK --help`

### Sample Output
```
$ transverseHFK -i 10 -X [10,5,8,6,3,7,2,4,9,1] -O [7,9,3,4,5,1,6,10,2,8]
    LL is null-homologous
    UR is NOT null-homologous
    D1[LL] is null-homologous
    D1[UR] is NOT null-homologous
```

## tHFK
To build and install the python library run
```
python setup.py install
```
or
```
make python-install
```

## Platform specific notes
Building on mac requires argp-standalone to be installed

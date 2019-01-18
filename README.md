# transverse-hfk-revision
This is a cleanup of the transverseHFK.c to compute the invariant defined in 
"Legendrian knots, transverse knots, and combinatorial Floer homology" by 
P. S. Ozsvath, Z. Szabo, and D. P. Thurston. The original code is from that 
linked in the paper "Transverse knots distinguished by Knot Floer Homology" 
by L. Ng, P. S. Ozsvath, and D. P. Thurston.

## Usage
Currently in order to change set the knot you must edit the file
src/TransverseHFK.c and set the lines
```
#define ArcIndex <grid size>
char Xs[ArcIndex] = {<Locations of the Xs seperated by commas>};
char Os[ArcIndex] = {<Locations of the Os seperated by commas>};
```
Go to the root directory of the repo and run
```
make
```
followed by
```
./transverseHFK
```
or the equivalent on your system. The program will calculate
whether or not x^-, x^+, delta_1(x^-), and delta_1(x^+) are null-homologous. 

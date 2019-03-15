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
To integrate tHFK with gridlink add the following lines to
the file `gridlink/gridlink.py` in the package

1. Add these lines to the list of imports

```
try:
    from tHFK import *
except:
    Tk_tHFK=None
```

2. Within the `__init__` method for the class Gridlink change this
```
if (TkHFK):
    invariantmenu = Menu(menubar, tearoff=0)
    invariantmenu.add_command(label='HFK^', command=self.HFKhat)
    menubar.add_cascade(label='Invariants', menu=invariantmenu)
```
to this
```
if (TkHFK):
    invariantmenu.add_command(label='HFK^', command=self.HFKhat)
if (Tk_tHFK):
    invariantmenu.add_command(label="transverseHFK", command=self.tHFK)
menubar.add_cascade(label='Invariants', menu=invariantmenu)
```

3. Add the following method to the gridlink class
```
def tHFK(self):
    if self.components > 1:
        showwarning('Knots only',
                    'Sorry, I can only compute transverseHFK invariants for knots')
        return
    Xlist, Olist = self.get_XOlists()
    t_hfk_object = Tk_tHFK([x+1 for x in Xlist], [o+1 for o in Olist], name=self.window.title())
```

Now after installing gridlink via setup.py and installing the tHFK python
library there will be a menu that allows you to run the transverseHFK 
invariants.

## Platform specific notes
Building on mac requires argp-standalone to be installed

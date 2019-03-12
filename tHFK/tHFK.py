from Tkinter import *
import ScrolledText
import _tHFK

class tHFK:
    """
    Class the contains the X,O code of a knot and
    provides methods for the transverseHFK invariants of
    the knot.

    Attributes
    ----------
    Xs : [int]
        A list of integers.
    Os : [int]
        A list of integers.

    Note: For the methods to work the Xs and Os must be
    permutations {1,...,N} with nonoverlapping values.

    Methods
    -------
    arc_index()
        Returns the size of the grid specified by the Xs and Os.
    x_plus()
        Returns an integer list corresponding to the x+ grid state.
    x_minus()
        Returns an integer list corresponding to the x- grid state.
    lambda_plus()
        Returns True if x+ is null-homologous. False otherwise.
    lambda_minus()
        Returns True if x- is null-homologous. False otherwise.
    d_lambda_plus()
        Returns True if d_1 x+ is null-homologous. False otherwise.
    d_lambda_minus()
        Returns True if d_1 x- is null-homologous. False otherwise.
    theta_n(n)
        Returns True if the lift of x+ to the n-fold cyclic branch cover
        is null-homologous. False otherwise.
    """

    def __init__(self, Xs, Os):
        """
        Parameters
        ----------
        Xs : [int]
           List of integers
        Os : [int]
           List of integers
        Note: For the methods to work the Xs and Os must be
        permutations {1,...,N} with nonoverlapping values.
        """
        self.Xs = Xs
        self.Os = Os

    def arc_index(self):
        """Returns the size of the grid specified by the Xs and Os."""
        return len(self.Xs)

    def x_plus(self):
        """Returns an integer list corresponding to the x+ grid state."""
        UR = [None]*self.arc_index()
        if self.Xs[self.arc_index()-1] == self.arc_index():
            UR[0] = 1
        else:
            UR[0] = self.Xs[self.arc_index()-1]+1
        for i in range(1,self.arc_index()):
            if self.Xs[i-1] == self.arc_index():
                UR[i] = 1
            else:
                UR[i] = self.Xs[i-1]+1
        return UR
                
    def x_minus(self):
        """Returns an integer list corresponding to the x- grid state."""
        return self.Xs
    
    def lambda_plus(self):
        """Returns True if x+ is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D0Q(self.x_plus(), self.Xs, self.Os)

    def lambda_minus(self):
        """Returns True if x- is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D0Q(self.x_minus(), self.Xs, self.Os)

    def d_lambda_plus(self):
        """Returns True if d_1 x+ is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D1Q(self.x_plus(), self.Xs, self.Os)

    def d_lambda_minus(self):
        """Returns True if d_1 x- is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D1Q(self.x_minus(), self.Xs, self.Os)

    def theta_n(self, n):
        """
        Returns True if the lift of x+ to the n-fold cyclic branch cover
        is null-homologous. False otherwise.
        """
        raise NotImplementedError
    
class Tk_tHFK(tHFK):
    """
    Tkinter window for use with the tHFK class methods

    Attributes
    ----------
    Xs : [int]
        A list of integers.
    Os : [int]
        A list of integers.
    name : str
        Name for the Tkinter window. Defaults to transverseHFK

    Note: For the methods to work the Xs and Os must be
    permutations {1,...,N} with nonoverlapping values.
    """
    def __init__(self,Xs,Os,name=None):
        """
        Initializes the Tkinter window.
        Parameters
        ----------
        Xs : [int]
            A list of integers.
        Os : [int]
            A list of integers.
        name : str
            Name for the Tkinter window. Defaults to transverseHFK

        Note: For the methods to work the Xs and Os must be
        permutations {1,...,N} with nonoverlapping values.
        """
        tHFK.__init__(self, Xs ,Os)
        self.window = Tk()
        if name:
            self.window.title(name)
        else:
            self.window.title("transverseHFK")

        self.l_plus_btn = Button(self.window, text=u"\u03BB^+", command=self.l_plus_btn_cmd)
        self.l_minus_btn = Button(self.window, text=u"\u03BB^-", command=self.l_minus_btn_cmd)
        self.d_plus_btn = Button(self.window, text=u"\u03B4_1 \u03BB^+", command=self.d_plus_btn_cmd)
        self.d_minus_btn = Button(self.window, text=u"\u03B4_1 \u03BB^-", command=self.d_minus_btn_cmd)
        self.theta_n_btn = Button(self.window, text=u"\u03B8_n", command=self.theta_n_btn_cmd)
        self.abort_btn = Button(self.window, text="Abort", command=self.abort_btn_cmd, state=DISABLED)
        self.n_lbl = Label(self.window,text="n=")
        self.n_var = StringVar()
        self.n_var.set("1")
        self.n_entry = Entry(self.window, width=3, text="n=")
        self.verbose_var = BooleanVar()
        self.verbose_var.set(False)
        self.verbosity_checkbox = Checkbutton(self.window, text="Verbose?", var=self.verbose_var)
        self.output_area = ScrolledText.ScrolledText(self.window,width=60,height=30)
        self.output_area.config(state=DISABLED)

        self.l_plus_btn.grid(column=0,row=0)
        self.l_minus_btn.grid(column=1,row=0)
        self.d_plus_btn.grid(column=2,row=0)
        self.d_minus_btn.grid(column=3,row=0)
        self.theta_n_btn.grid(column=4,row=0)
        self.abort_btn.grid(column=2,row=3)
        self.n_lbl.grid(column=5,row=0)
        self.n_entry.grid(column=6,row=0)
        self.verbosity_checkbox.grid(column=0,row=1)
        self.output_area.grid(column=0,row=2, columnspan=7)

        self.window.mainloop()

    def _writeln_output(self,s):
        """
        Writes the string s to the output area.
        Parameters
        ----------
        s : str
        """
        self.output_area.config(state=NORMAL)
        self.output_area.insert(END,'\n')
        self.output_area.insert(END,s)
        self.output_area.config(state=DISABLED)

    def _with_abort(self,f,*args):
        """
        Calls f with arguments args such that the abort button is active.
        
        Parameters
        ----------
        f: fun
            The function that should have abort functionality
        args:
            List containing arguments for the function f
        """
        self.abort_btn.config(state=NORMAL)
        result = f(*args)
        self.abort_btn.config(state=DISABLED)
        return result
        
    def l_plus_btn_cmd(self):
        """Calls lambda_plus and prints the result to the output area."""
        if self._with_abort(self.lambda_plus):
            self._writeln_output(u"\u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03BB^+ is NOT null-homologous")

    def l_minus_btn_cmd(self):
        """Calls lambda_minus and prints the result to the output area."""
        if self._with_abort(self.lambda_minus):
            self._writeln_output(u"\u03BB^- is null-homologous")
        else:
            self._writeln_output(u"\u03BB^- is NOT null-homologous")

    def d_plus_btn_cmd(self):
        """Calls d_lambda_plus and prints the result to the output area."""
        if self._with_abort(self.d_lambda_plus):
            self._writeln_output(u"\u03B4_1 \u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03B4_1 \u03BB^+ is NOT null-homologous")

    def d_minus_btn_cmd(self):
        """Calls d_lambda_minus and prints the result to the output area."""
        if self._with_abort(self.d_lambda_minus):
            self._writeln_output(u"\u03B4_1 \u03BB^- is null-homologous")
        else:
            self._writeln_output(u"\u03B4_1 \u03BB^- is NOT null-homologous")

    def theta_n_btn_cmd(self):
        """
        Calls theta_n(n) with n supplied from the entry widget and
        prints the result to the output area.
        """
        try:
            n = int(self.n_entry.get())
        except ValueError:
            self._writeln_output("Error: n must be an integer")
            return
        
        if self._with_abort(self.theta_n,n):
            self._writeln_output(u"\u03B8_n is null-homologous")
        else:
            self._writeln_output(u"\u03B8_n is NOT null-homologous")

    def abort_btn_cmd(self):
        """Halts the currently running execution"""
        raise NotImplementedError

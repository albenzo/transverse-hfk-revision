from Tkinter import *
import ScrolledText
from sys import stdout
import multiprocessing as mp
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

    def __init__(self, Xs, Os, out_stream=stdout, verbosity=0):
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
        self.out_stream = out_stream
        self.verbosity = verbosity

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
        return _tHFK.null_homologous_D0Q(self.x_plus(), self.Xs, self.Os, self.out_stream, self.verbosity)

    def lambda_minus(self):
        """Returns True if x- is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D0Q(self.x_minus(), self.Xs, self.Os, self.out_stream, self.verbosity)

    def d_lambda_plus(self):
        """Returns True if d_1 x+ is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D1Q(self.x_plus(), self.Xs, self.Os, self.out_stream, self.verbosity)

    def d_lambda_minus(self):
        """Returns True if d_1 x- is null-homologous. False otherwise."""
        return _tHFK.null_homologous_D1Q(self.x_minus(), self.Xs, self.Os, self.out_stream, self.verbosity)

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
        tHFK.__init__(self, Xs ,Os, self, 0)
        self.window = Tk()
        if name:
            self.name = name
            self.window.title(name)
        else:
            self.name = "transverseHFK"
            self.window.title("transverseHFK")

        self._process_list = []
        self._write_queue = mp.Queue()
        
        self.output_area = ScrolledText.ScrolledText(self.window,width=60,height=30)
        self.output_area.config(state=DISABLED)
        self.l_plus_btn = Button(self.window, text=u"\u03BB^+", command=self._with_process(self.l_plus_btn_cmd))
        self.l_minus_btn = Button(self.window, text=u"\u03BB^-", command=self._with_process(self.l_minus_btn_cmd))
        self.d_plus_btn = Button(self.window, text=u"\u03B4_1 \u03BB^+", command=self._with_process(self.d_plus_btn_cmd))
        self.d_minus_btn = Button(self.window, text=u"\u03B4_1 \u03BB^-", command=self._with_process(self.d_minus_btn_cmd))
        self.n_lbl = Label(self.window,text="n=")
        self.n_entry = Spinbox(self.window, width=3, text="n=", from_=1, to=100)
        self.theta_n_btn = Button(self.window, text=u"\u03B8_n", command=self._with_process(self.theta_n_btn_cmd), state=DISABLED)
        self.abort_btn = Button(self.window, text="Abort", command=self.abort_btn_cmd)
        self.verbosity_var = StringVar()
        self.verbosity_var.trace("w", self._sync_verbosity)
        self.verbosity_var.set('silent')
        self.verbosity_list = OptionMenu(self.window, self.verbosity_var, 'silent', 'quiet', 'verbose')
                
        self.l_plus_btn.grid(column=0,row=0)
        self.l_minus_btn.grid(column=1,row=0)
        self.d_plus_btn.grid(column=2,row=0)
        self.d_minus_btn.grid(column=3,row=0)
        self.theta_n_btn.grid(column=4,row=0)
        self.abort_btn.grid(column=2,row=3)
        self.n_lbl.grid(column=5,row=0)
        self.n_entry.grid(column=6,row=0)
        self.verbosity_list.grid(column=0,row=1)
        self.output_area.grid(column=0,row=2, columnspan=7)

        self._callback_id = self.window.after(0,self._queue_check)
        self.window.protocol("WM_DELETE_WINDOW", self._clean_and_destroy)
        self.window.mainloop()

    def _clean_and_destroy(self):
        """Ends _write_queue polling and destroys the window"""
        self.window.after_cancel(self._callback_id)
        self.window.destroy()

    def _sync_verbosity(self, *args):
        self.verbosity = ['silent','quiet','verbose'].index(self.verbosity_var.get())
        
    def write(self,s):
        """
        Sends the string s to _write_queue where it will be processed
        and written to output_area by the event loop.
        
        Parameters
        ----------
        s : str
        """
        self._write_queue.put(s)

    def _write_output_area(self,s):
        """
        Writes the string s to the output_area.

        Parameters
        ----------
        s : str
        """
        self.output_area.config(state=NORMAL)
        self.output_area.insert(END,s)
        self.output_area.config(state=DISABLED)


    def _with_process(self,f):
        """
        Takes in a no parameter function f and returns a closure
        that runs f in a new daemon.

        Parameters
        ----------
        f : fun: () -> ()
        """
        def p_f():
            p = mp.Process(target=f)
            self._process_list.append(p)
            p.daemon = True
            p.start()
        return p_f

    def _queue_check(self):
        """
        Polls _write_queue without blocking and schedules any strings
        to be written.
        """
        while True:
            try:
                s = self._write_queue.get_nowait()
            except mp.queues.Empty:
                break
            else:
                #self._write_output_area(s)
                self.window.after_idle(self._write_output_area, s)
        self._callback_id = self.window.after(100, self._queue_check)
    
    def l_plus_btn_cmd(self):
        """Calls lambda_plus and prints the result to the output area."""
        try:
            if self.lambda_plus():
                self.write(u"\u03BB^+ is null-homologous\n")
            else:
                self.write(u"\u03BB^+ is NOT null-homologous\n")
        except Exception as e:
            self.write(str(e))

    def l_minus_btn_cmd(self):
        """Calls lambda_minus and prints the result to the output area."""
        try:
            if self.lambda_minus():
                self.write(u"\u03BB^- is null-homologous\n")
            else:
                self.write(u"\u03BB^- is NOT null-homologous\n")
        except Exception as e:
            self.write(str(e))

    def d_plus_btn_cmd(self):
        """Calls d_lambda_plus and prints the result to the output area."""
        try:
            if self.d_lambda_plus():
                self.write(u"\u03B4_1 \u03BB^+ is null-homologous\n")
            else:
                self.write(u"\u03B4_1 \u03BB^+ is NOT null-homologous\n")
        except Exception as e:
            self.write(str(e))

    def d_minus_btn_cmd(self):
        """Calls d_lambda_minus and prints the result to the output area."""
        try:
            if self.d_lambda_minus():
                self.write(u"\u03B4_1 \u03BB^- is null-homologous\n")
            else:
                self.write(u"\u03B4_1 \u03BB^- is NOT null-homologous\n")
        except Exception as e:
            self.write(str(e))

    def theta_n_btn_cmd(self):
        """
        Calls theta_n(n) with n supplied from the entry widget and
        prints the result to the output area.
        """
        try:
            n = int(self.n_entry.get())
        except ValueError:
            self.write("Error: n must be an integer\n")
            return
        
        try:
            if self.theta_n(n):
                self.write(u"\u03B8_" + str(n) + " is null-homologous\n")
            else:
                self.write(u"\u03B8_" + str(n) + " is NOT null-homologous\n")
        except Exception as e:
            self.write(str(e))

    def abort_btn_cmd(self):
        """
        Sends SIGTERM to any processes spawned by the window. Then
        creates a new _write_queue.
        """
        for p in self._process_list:
            if p.is_alive():
                p.terminate()
        self._process_list = []
        self._write_queue = mp.Queue()

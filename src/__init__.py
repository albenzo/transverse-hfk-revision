import tHFK._tHFK as _tHFK
from Tkinter import *
import ScrolledText

__all__ = ['tHFK', 'Tk_tHFK']

class tHFK:
    """
    """

    def __init__(self, Xs, Os, name=None):
        """
        """
        self.Xs = Xs
        self.Os = Os
        self.name = name

    def x_plus(self):
        UR = [None]*len(self.Xs)
        for i in range(len(self.Xs)):
            if(self.Xs[len(self.Xs)-1] == len(self.Xs)):
                UR[0] = 1
            else:
                UR[i] = (self.Xs[i-1]+1) % len(self.Xs)
        return [i+1 for i in UR]
    
    def x_minus(self):
        return self.Xs
    
    def lambda_plus(self):
        print self.Xs
        print self.Os
        print self.x_plus()
        print self.x_minus()
        return _tHFK.null_homologous_D0Q(self.x_plus(), self.Xs, self.Os)

    def lambda_minus(self):
        return _tHFK.null_homologous_D0Q(self.x_minus(), self.Xs, self.Os)

    def d_lambda_plus(self):
        return _tHFK.null_homologous_D1Q(self.x_plus(), self.Xs, self.Os)

    def d_lambda_minus(self):
        return _tHFK.null_homologous_D1Q(self.x_minus(), self.Xs, self.Os)

    def theta_n(self, n):
        raise NotImplementedError
    
class Tk_tHFK(tHFK):
    """
    """
    def __init__(self,Xs,Os,name=None):
        tHFK.__init__(self, Xs ,Os, name)
        self.window = Tk()
        if self.name:
            self.window.title(self.name)
        else:
            self.window.title("transverseHFK")

        self.l_plus_btn = Button(self.window, text=u"\u03BB^+", command=self.l_plus_btn_cmd)
        self.l_minus_btn = Button(self.window, text=u"\u03BB^-", command=self.l_minus_btn_cmd)
        self.dl_plus_btn = Button(self.window, text=u"\u03B4_1 \u03BB^+", command=self.dl_plus_btn_cmd)
        self.dl_minus_btn = Button(self.window, text=u"\u03B4_1 \u03BB^-", command=self.dl_minus_btn_cmd)
        self.theta_n_btn = Button(self.window, text=u"\u03B8_n", command=self.theta_n_btn_cmd)
        self.abort_btn = Button(self.window, text="Abort", command=self.abort_btn_cmd)
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
        self.dl_plus_btn.grid(column=2,row=0)
        self.dl_minus_btn.grid(column=3,row=0)
        self.theta_n_btn.grid(column=4,row=0)
        self.abort_btn.grid(column=2,row=3)
        self.n_lbl.grid(column=5,row=0)
        self.n_entry.grid(column=6,row=0)
        self.verbosity_checkbox.grid(column=0,row=1)
        self.output_area.grid(column=0,row=2, columnspan=7)

        self.window.mainloop()

    def _writeln_output(self,s):
        self.output_area.config(state=NORMAL)
        self.output_area.insert(END,'\n')
        self.output_area.insert(END,s)
        self.output_area.config(state=DISABLED)
        
    def l_plus_btn_cmd(self):
        if self.lambda_plus():
            self._writeln_output(u"\u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03BB^+ is NOT null-homologous")

    def l_minus_btn_cmd(self):
        if self.lambda_minus():
            self._writeln_output(u"\u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03BB^+ is NOT null-homologous")

    def dl_plus_btn_cmd(self):
        if self.d_lambda_plus():
            self._writeln_output(u"\u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03BB^+ is NOT null-homologous")

    def dl_minus_btn_cmd(self):
        if self.d_lambda_minus():
            self._writeln_output(u"\u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03BB^+ is NOT null-homologous")

    def theta_n_btn_cmd(self):
        try:
            n = int(self.n_entry.get())
        except ValueError:
            self._writeln_output("Error: n must be an integer")
            return
        
        if self.theta_n(n):
            self._writeln_output(u"\u03BB^+ is null-homologous")
        else:
            self._writeln_output(u"\u03BB^+ is NOT null-homologous")

    def abort_btn_cmd(self):
        raise NotImplementedError

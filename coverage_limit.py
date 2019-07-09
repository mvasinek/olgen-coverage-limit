from tkinter import *
from binom import *

class Window(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master
        self.result_coverage = StringVar()
        self.result_variant_reads = StringVar()

        self.err_label = None

        self.params_coverage = (("vaf", "Variant allele frequency (VAF)",""),
                                ("err", "Sequencing error","0.01"),
                                ("err_lim", "Limit of cumulative probability for errors\n(i.e. 1 - probability of false positive)","0.001"),
                                ("det_lim", "Limit of cumulative probability for number of variant reads\n(i.e. probability of true positive)","0.999"),
                                ("min_var", "Minimal number of required variant reads\n(optional parameter)","")
        )

        self.params_variant_reads = (("cov", "Coverage depth",""),
                                     ("err", "Sequencing error",""),
                                     ("err_lim", "Limit of cumulative probability error","0.95"),
        )
        
        self.entries_coverage = {}
        self.entries_variant_reads = {}
        
        self.initWindow()

    def initWindow(self):
        self.master.title("OLGEN Coverage Limit")
        self.pack(fill=BOTH, expand=1)

        parent_panel = PanedWindow(self, orient=HORIZONTAL)
        parent_panel.pack(fill=BOTH, expand=0)
        
        coverage_frame = LabelFrame(parent_panel, text="Coverage")
        coverage_frame.pack(fill=BOTH, expand=1)
        #parent_panel.add(coverage_frame)

        coverage_item_frame = Frame(coverage_frame,pady=10)
        coverage_item_frame.pack()
        """
        COVERAGE
        """
        r = 0
        for param in self.params_coverage:
            item_label = Label(coverage_item_frame, text=param[1])
            item_label.grid(row=r,column=0)
            item_entry = Entry(coverage_item_frame)
            item_entry.grid(row=r,column=1)
            item_entry.bind('<Return>', self.computeCoverage)
            item_entry.insert(END, param[2])
            if r == 0:
                item_entry.focus()

            self.entries_coverage[param[0]] = item_entry

            r += 1

        coverage_start_button = Button(coverage_frame, text="Compute coverage", command=self.computeCoverage)
        coverage_start_button.pack(side=BOTTOM)

        coverage_result_label = Label(coverage_frame, textvariable=self.result_coverage, font=('Helvetica', 10, 'bold'))
        coverage_result_label.pack(side=BOTTOM)

        self.err_label = coverage_result_label

        variant_result_label = Label(coverage_frame, textvariable=self.result_variant_reads, font=('Helvetica', 10, 'bold'))
        variant_result_label.pack(side=BOTTOM)

    def computeCoverage(self, component=None):
        self.result_coverage.set("")
        self.err_label.configure(fg="black")
        self.result_variant_reads.set("")

        vaf = 0
        err = 0
        err_lim = 0
        det_lim = 0
        min_var = None

        try:
            vaf = float(self.entries_coverage['vaf'].get())
        except:
            self.result_coverage.set("Cannot convert variant allele frequency\nto float.")
            self.err_label.configure(fg='red')
            return
        
        if vaf <= 0.0 or vaf > 1.0:
            self.result_coverage.set("VAF must be in (0;1).")
            self.err_label.configure(fg='red')
            return

        try:
            err = float(self.entries_coverage['err'].get())
        except:
            self.result_coverage.set("Cannot convert sequencing error to float.")
            self.err_label.configure(fg='red')
            return

        if err < 0.0 or err > 1.0:
            self.result_coverage.set("Sequencing error must be in <0;1).")
            self.err_label.configure(fg='red')
            return
        
        if vaf <= err:
            self.result_coverage.set("VAF must exceed sequencing error.")
            self.err_label.configure(fg='red')
            return

        try:
            err_lim = float(self.entries_coverage['err_lim'].get())
        except:
            self.result_coverage.set("Cannot convert limit of cumulative\nprobability error to float.")
            self.err_label.configure(fg='red')
            return
        
        if err_lim < 0.0 or err_lim > 1.0:
            self.result_coverage.set("Cumulative probability limit for\ncoverage must be in <0;1>.")
            self.err_label.configure(fg='red')
            return

        try:
            det_lim = float(self.entries_coverage['det_lim'].get())
        except:
            self.result_coverage.set("Cannot convert limit of cumulative\nprobability coverage to float.")
            self.err_label.configure(fg='red')
            return        

        if err_lim == 1.0 and det_lim == 1.0:
            self.result_coverage.set("Having both limits set to 1.0 would\nnever converge.")
            self.err_label.configure(fg='red')
            return

        if det_lim < 0.0 or det_lim > 1.0:
            self.result_coverage.set("Cumulative probability limit for sequencing\nerror must be in <0;1>.")
            self.err_label.configure(fg='red')
            return

        try:
            min_var_val = self.entries_coverage['min_var'].get()
            if min_var_val != "":
                min_var = int(min_var_val)
        except:
            self.result_coverage.set("Cannot convert minimal number of variant\n reads to integer.")
            self.err_label.configure(fg='red')
            return        

        if min_var == None:
            cov,variants = relaxTPFP(INITIAL_COVERAGE, vaf, err, 1-err_lim, det_lim)
        else:
            cov,variants = relaxTPFPMin(INITIAL_COVERAGE, vaf, err, 1-err_lim, det_lim, min_var)

        self.result_coverage.set("Recommended coverage: %d" % (cov))
        self.result_variant_reads.set("Minimal number of variant reads: %d" % (variants))

root = Tk()
app = Window(root)
root.mainloop()

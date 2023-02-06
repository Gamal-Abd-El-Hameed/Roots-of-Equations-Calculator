import datetime
import time
import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext

import matplotlib.pyplot as plt
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from sympy import *
from Bracketing.Bisection import Bisection
from Bracketing.Bracketing import Bracketing
from Bracketing.FalsePosition import FalsePosition
from FixedPoint.FixedPoint import FixedPointIteration
from Secant.Secant import Secant
from NewtonRaphson.NewtonRaphson import NewtonRaphson

from GaussJaccobi.GaussSeidil import GaussSeidil
from GaussJaccobi.Jacobi import Jacobi
from Gauss.GaussEliminator import GaussEliminator
from GJ.FloatChopper import FloatChopper
from GJ.FloatRounder import FloatRounder
from GJ.GaussJordan import Gauss_Jordan

from inputHandler import *
from LU.LUController import *


class App(tk.Tk):
    def __init__(self):
        super().__init__()

        self.bracketing_guess_entry = None
        self.bracketing_delta_entry = None
        self.l10 = None
        self.l9 = None
        self.bracketing_options = None
        self.open_method_2_entry = None
        self.open_method_1_entry = None
        self.l8 = None
        self.l7 = None
        self.equations_input_1 = None
        self.l0_1 = None
        self.l6 = None
        self.l5 = None
        self.bracketing_entry_Xl = None
        self.bracketing_entry_Xu = None
        self.initial_guess_entry = None
        self.l4 = None
        self.l3 = None
        self.l2 = None
        self.l1 = None
        self.l0 = None
        self.epsilon_entry = None
        self.iterations_used = None
        self.lu_options = None
        self.digits_used = None
        self.equations_input = None
        self.method_menu = None
        self.submit_button = None

        self.title('Numerical Project')
        self.geometry("1300x410")

        self.number_of_iterations = tk.StringVar()
        self.number_of_iterations.set("50")

        self.epsilon = tk.StringVar()
        self.epsilon.set("10**-5")

        self.lu_config = tk.StringVar()
        self.lu_config.set("Doolittle")
        
        self.bracketing_config = tk.StringVar()
        self.bracketing_config.set("Xl & Xu")

        self.x_lower = tk.StringVar()
        self.x_lower.set("")

        self.x_upper = tk.StringVar()
        self.x_upper.set("")
        
        self.bracketing_delta = tk.StringVar()
        self.bracketing_delta.set("")
        
        self.bracketing_guess = tk.StringVar()
        self.bracketing_guess.set("")
        
        self.first_guess = tk.StringVar()
        self.first_guess.set("")

        self.second_guess = tk.StringVar()
        self.second_guess.set("")

        self.number_of_digits = tk.StringVar()
        self.number_of_digits.set("7")

        self.methods = ["Gauss", "Gauss-Jordan", "LU", "Gauss-Seidel", "Jacobi", "Bisection", "False-Position", "Fixed point", "Newton-Raphson", "Secant"]
        self.operation_on_numbers = ["Rounding", "Chopping"]

        self.round_or_chop = tk.StringVar()
        self.round_or_chop.set("Operation on numbers")

        self.method_selected = tk.StringVar()
        self.method_selected.set("Select the method")
        # Define the grid
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=1)
        self.columnconfigure(4, weight=1)
        self.columnconfigure(5, weight=1)
        self.columnconfigure(6, weight=1)
        self.columnconfigure(7, weight=1)
        self.columnconfigure(8, weight=1)
        self.columnconfigure(9, weight=1)
        self.columnconfigure(10, weight=1)
        self.columnconfigure(11, weight=1)
        self.columnconfigure(12, weight=1)
        self.columnconfigure(13, weight=1)
        self.columnconfigure(14, weight=1)
        self.columnconfigure(15, weight=1)
        self.columnconfigure(16, weight=1)
        self.columnconfigure(17, weight=1)
        self.columnconfigure(18, weight=1)
        self.columnconfigure(19, weight=1)
        self.columnconfigure(20, weight=1)

        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)
        self.rowconfigure(4, weight=1)
        self.rowconfigure(5, weight=1)
        self.rowconfigure(6, weight=1)
        self.rowconfigure(7, weight=1)
        self.rowconfigure(8, weight=1)
        self.rowconfigure(9, weight=1)
        self.rowconfigure(10, weight=1)
        self.rowconfigure(11, weight=1)
        self.rowconfigure(12, weight=1)
        self.rowconfigure(13, weight=1)
        self.rowconfigure(14, weight=1)
        self.rowconfigure(15, weight=1)
        self.rowconfigure(16, weight=1)
        self.rowconfigure(17, weight=1)
        self.rowconfigure(18, weight=1)
        # Creating The widgets
        self.create_widgets()

    def create_widgets(self):
        self.method_menu = ttk.OptionMenu(self, self.method_selected, self.method_selected.get(), *self.methods, command=self.toggle)
        self.method_menu.grid(column=0, row=0)

        self.method_menu = ttk.OptionMenu(self, self.round_or_chop, self.round_or_chop.get(), *self.operation_on_numbers)
        self.method_menu.grid(column=1, row=0)

        self.lu_options = ttk.OptionMenu(self, self.lu_config, self.lu_config.get(), *["Doolittle", "Crout"])
        self.bracketing_options = ttk.OptionMenu(self, self.bracketing_config, self.bracketing_config.get(), *["Xl & Xu"], command=self.bracketing_changed)
        
        self.l1 = ttk.Label(self, text="Digits:")
        self.l1.grid(column=19, row=0)
        self.digits_used = ttk.Entry(self, textvariable=self.number_of_digits, justify='center')
        self.digits_used.grid(column=20, row=0)

        self.l2 = ttk.Label(self, text="Iterations:")
        self.iterations_used = ttk.Entry(self, textvariable=self.number_of_iterations, justify='center')

        self.l3 = ttk.Label(self, text="Epsilon:")
        self.epsilon_entry = ttk.Entry(self, textvariable=self.epsilon, justify='center')

        self.l4 = ttk.Label(self, text="Initial guess", font=10)
        self.initial_guess_entry = scrolledtext.ScrolledText(self, wrap="none", width=25, height=10)

        self.l5 = ttk.Label(self, text="X Lower:")
        self.bracketing_entry_Xl = ttk.Entry(self, textvariable=self.x_lower, justify='center')

        self.l6 = ttk.Label(self, text="X Upper:")
        self.bracketing_entry_Xu = ttk.Entry(self, textvariable=self.x_upper, justify='center')

        self.l7 = ttk.Label(self, text="Initial guess (X0):")
        self.open_method_1_entry = ttk.Entry(self, textvariable=self.first_guess, justify='center')

        self.l8 = ttk.Label(self, text="Initial guess (X1):")
        self.open_method_2_entry = ttk.Entry(self, textvariable=self.second_guess, justify='center')

        self.l9 = ttk.Label(self, text="Initial guess:")
        self.bracketing_guess_entry = ttk.Entry(self, textvariable=self.bracketing_guess, justify='center')

        self.l10 = ttk.Label(self, text="Delta:")
        self.bracketing_delta_entry = ttk.Entry(self, textvariable=self.bracketing_delta, justify='center')

        self.l0 = ttk.Label(self, text="Equations", font=10)
        self.l0.grid(column=10, row=4)

        self.l0_1 = ttk.Label(self, text="Enter an equation", font=10)

        self.equations_input = scrolledtext.ScrolledText(self, wrap="none", width=60, height=15)
        self.equations_input.grid(column=10, row=5)

        self.equations_input_1 = scrolledtext.ScrolledText(self, wrap="none", width=60, height=1)

        self.submit_button = ttk.Button(self, text="Solve", command=self.solve_equations)
        self.submit_button.grid(column=10, row=7)

    def toggle(self, *args):
        self.l0.grid(column=10, row=4)
        self.equations_input.grid(column=10, row=5)

        self.first_guess.set("")
        self.second_guess.set("")
        self.bracketing_guess.set("")
        self.bracketing_delta.set("")
        self.x_upper.set("")
        self.x_lower.set("")
        self.bracketing_config.set("Xl & Xu")
        self.lu_config.set("Doolittle")

        if self.number_of_digits.get() == "":
            self.number_of_digits.set("7")

        if self.epsilon.get() == "":
            self.epsilon.set("10**-5")

        if self.number_of_iterations.get() == "":
            self.number_of_iterations = "50"

        self.l0_1.grid_forget()
        self.l2.grid_forget()
        self.l3.grid_forget()
        self.l4.grid_forget()
        self.l5.grid_forget()
        self.l6.grid_forget()
        self.l7.grid_forget()
        self.l8.grid_forget()
        self.l9.grid_forget()
        self.l10.grid_forget()
       
        self.iterations_used.grid_forget()
        self.epsilon_entry.grid_forget()
        self.initial_guess_entry.grid_forget()
        self.lu_options.grid_forget()
        self.bracketing_options.grid_forget()
        self.bracketing_entry_Xl.grid_forget()
        self.bracketing_entry_Xu.grid_forget()
        self.equations_input_1.grid_forget()
        self.open_method_1_entry.grid_forget()
        self.open_method_2_entry.grid_forget()
        self.bracketing_delta_entry.grid_forget()
        self.bracketing_guess_entry.grid_forget()

        root_finding_methods = ["Bisection", "False-Position", "Fixed point", "Newton-Raphson", "Secant"]
        try:
            if self.method_selected.get() in root_finding_methods:
                self.l0.grid_forget()
                self.l0_1.grid(column=10, row=4)
                self.equations_input.grid_forget()
                self.equations_input_1.grid(column=10, row=5)
                self.l2.grid(column=19, row=1)
                self.iterations_used.grid(column=20, row=1)
                self.l3.grid(column=19, row=2)
                self.epsilon_entry.grid(column=20, row=2)
        except:
            print('', end='')

        if self.method_selected.get() == "LU":
            self.lu_options.grid(column=2, row=0)

        elif self.method_selected.get() == "Gauss-Seidel" or self.method_selected.get() == "Jacobi":
            self.l2.grid(column=19, row=1)
            self.iterations_used.grid(column=20, row=1)
            self.l3.grid(column=19, row=2)
            self.epsilon_entry.grid(column=20, row=2)
            self.l4.grid(column=19, row=4)
            self.initial_guess_entry.grid(column=19, row=5)

        elif self.method_selected.get() == "Bisection" or self.method_selected.get() == "False-Position":
            self.bracketing_options.grid(column=2, row=0)
            self.l5.grid(column=17, row=0)
            self.bracketing_entry_Xl.grid(column=18, row=0)
            self.l6.grid(column=17, row=1)
            self.bracketing_entry_Xu.grid(column=18, row=1)

        elif self.method_selected.get() == "Fixed point" or self.method_selected.get() == "Newton-Raphson":
            self.l7.grid(column=17, row=0)
            self.open_method_1_entry.grid(column=18, row=0)
            
        elif self.method_selected.get() == "Secant":
            self.l7.grid(column=17, row=0)
            self.open_method_1_entry.grid(column=18, row=0)
            self.l8.grid(column=17, row=1)
            self.open_method_2_entry.grid(column=18, row=1)

    def bracketing_changed(self, *args):
        if self.bracketing_config.get() == "Xl & Xu":
            self.bracketing_guess.set("")
            self.bracketing_delta.set("")
            self.l9.grid_forget()
            self.l10.grid_forget()
            self.bracketing_guess_entry.grid_forget()
            self.bracketing_delta_entry.grid_forget()

            self.l5.grid(column=17, row=0)
            self.bracketing_entry_Xl.grid(column=18, row=0)
            self.l6.grid(column=17, row=1)
            self.bracketing_entry_Xu.grid(column=18, row=1)
        elif self.bracketing_config.get() == "Initial guess":
            self.x_upper.set("")
            self.x_lower.set("")
            self.l5.grid_forget()
            self.bracketing_entry_Xl.grid_forget()
            self.l6.grid_forget()
            self.bracketing_entry_Xu.grid_forget()

            self.l9.grid(column=17, row=0)
            self.bracketing_guess_entry.grid(column=18, row=0)
            self.l10.grid(column=17, row=1)
            self.bracketing_delta_entry.grid(column=18, row=1)

    def isfloat(self, num):
        try:
            float(num)
            return True
        except ValueError:
            return False

    def solve_equations(self):

        def close_error():
            error_window.destroy()

        def close_answer():
            answer_window.destroy()

        if self.epsilon.get() == "":
            self.epsilon.set("10**-5")

        eps = 10 ** -5
        initial_array = []

        if self.method_selected.get() == "Gauss-Seidel" or self.method_selected.get() == "Jacobi":
            try:
                eps = eval(self.epsilon.get())
            except:
                eps = 10 ** -5
                self.epsilon.set("10**-5")

        root_finding_methods = ["Bisection", "False-Position", "Fixed point", "Newton-Raphson", "Secant"]

        if self.method_selected.get() in root_finding_methods:
            equations = self.equations_input_1.get('1.0', 'end')
        else:
            equations = self.equations_input.get('1.0', 'end')

        initial = self.initial_guess_entry.get('1.0', 'end')

        if equations == "\n":
            error_window = tk.Toplevel(self)
            error_window.geometry("400x100")
            error_window.title("Error!")
            l1 = ttk.Label(error_window, text=f"Error!", font=14)
            l1.pack()
            l2 = ttk.Label(error_window, text=f"Enter The equations", font=12)
            l2.pack()
            b1 = ttk.Button(error_window, text="OK", command=close_error)
            b1.pack()
            return

        if self.method_selected.get() == "Bisection" or self.method_selected.get() == "False-Position":
            if ((self.x_lower.get() == "" or self.x_upper.get() == "") and self.bracketing_config.get() == "Xl & Xu") or ((self.bracketing_guess.get() == "" or self.bracketing_delta.get() == "") and self.bracketing_config.get() == "Initial guess"):
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l1.pack()
                l2 = ttk.Label(error_window, text=f"Check The bracketing input", font=12)
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return
            if ((not self.isfloat(self.x_lower.get()) or not self.isfloat(self.x_lower.get())) and self.bracketing_config.get() == "Xl & Xu") or ((not self.isfloat(self.bracketing_guess.get()) or not self.isfloat(self.bracketing_delta.get())) and self.bracketing_config.get() == "Initial guess"):
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l1.pack()
                l2 = ttk.Label(error_window, text=f"Check The bracketing input", font=12)
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return
        if self.method_selected.get() == "Fixed point" or self.method_selected.get() == "Newton-Raphson":
            if self.first_guess.get() == "" or (not self.isfloat(self.first_guess.get())):
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l1.pack()
                l2 = ttk.Label(error_window, text=f"Check The guessing input", font=12)
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        if self.method_selected.get() == "Secant":
            if self.first_guess.get() == "" or self.second_guess.get() == "" or (not self.isfloat(self.first_guess.get())) or (not self.isfloat(self.second_guess.get())):
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l1.pack()
                l2 = ttk.Label(error_window, text=f"Check The guessing input", font=12)
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        if not self.number_of_iterations.get().isnumeric():
            self.number_of_iterations.set("50")

        if not self.number_of_digits.get().isnumeric():
            self.number_of_digits.set("7")

        if self.method_selected.get() == "Gauss-Seidel" or self.method_selected.get() == "Jacobi":
            if initial == "\n":
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l1.pack()
                l2 = ttk.Label(error_window, text=f"Enter The initial values", font=12)
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

            initial_array = str.split(initial, "\n")
            i = len(initial_array) - 1
            while i >= 0:
                if initial_array[i] == "":
                    del initial_array[i]
                i = i - 1

            for k in range(len(initial_array)):
                try:
                    initial_array[k] = float(initial_array[k])
                except:
                    error_window = tk.Toplevel(self)
                    error_window.geometry("400x100")
                    error_window.title("Error!")
                    l1 = ttk.Label(error_window, text=f"Error!", font=14)
                    l2 = ttk.Label(error_window, text=f"Enter initial values correctly", font=12)
                    l1.pack()
                    l2.pack()
                    b1 = ttk.Button(error_window, text="OK", command=close_error)
                    b1.pack()
                    return

        if self.method_selected.get() == "LU":
            if self.lu_config.get() == "None":
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"Choose from a LU-Method", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        if self.method_selected.get() == "Select the method" or self.round_or_chop.get() == "Operation on numbers":
            error_window = tk.Toplevel(self)
            error_window.geometry("400x100")
            error_window.title("Error!")
            l1 = ttk.Label(error_window, text=f"Error!", font=14)
            l2 = ttk.Label(error_window, text=f"Choose from option-menu", font=12)
            l1.pack()
            l2.pack()
            b1 = ttk.Button(error_window, text="OK", command=close_error)
            b1.pack()
            return

        answer_window = tk.Toplevel(self)
        answer_window.geometry("600x400")
        answer_window.title("Solution")

        answer_window.columnconfigure(0, weight=1)
        answer_window.columnconfigure(1, weight=1)
        answer_window.columnconfigure(2, weight=1)
        answer_window.columnconfigure(3, weight=1)
        answer_window.columnconfigure(4, weight=1)
        answer_window.columnconfigure(5, weight=1)
        answer_window.columnconfigure(6, weight=1)

        answer_window.rowconfigure(0, weight=1)
        answer_window.rowconfigure(1, weight=1)
        answer_window.rowconfigure(2, weight=1)
        answer_window.rowconfigure(3, weight=1)
        answer_window.rowconfigure(4, weight=1)
        answer_window.rowconfigure(5, weight=1)
        answer_window.rowconfigure(6, weight=1)

        l1 = ttk.Label(answer_window, text=f"Solved By {self.method_selected.get()}", font=14)
        l1.grid(column=3, row=0)
        coefficients = []
        variables = []
        constants = []
        equation_str = ""
        if not (self.method_selected.get() in root_finding_methods):
            equations_array = str.split(equations, "\n")
            try:
                input_handler = InputHandler(equations_array)
                coefficients, variables, constants = input_handler.getAug()
            except:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"Check the equations", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return
        else:
            eqn_one = str(equations).split('\n')
            equation_str = str(eqn_one[0]).replace('\n', '').replace(' ', '').replace('^', '**')

        time_before = time.perf_counter_ns()

        answer = []
        the_plot = None

        if self.round_or_chop.get() == "Rounding":
            float_converter = FloatRounder(int(self.number_of_digits.get()))
        elif self.round_or_chop.get() == "Chopping":
            float_converter = FloatChopper(int(self.number_of_digits.get()))
        else:
            float_converter = None

        if self.method_selected.get() == "Gauss":
            try:
                gauss = GaussEliminator(coefficients, constants, float_converter)
                answer = gauss.solve()
            except:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"Check the input", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "Gauss-Jordan":
            try:
                gauss_jordan = Gauss_Jordan(coefficients, constants, float_converter)
                answer = gauss_jordan.solve()
            except:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"Check the input", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "LU":
            try:
                lu = LUController(self.lu_config.get(), coefficients, constants, float_converter)
                answer = lu.solve()
            except:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"Check the input", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "Gauss-Seidel":
            try:
                gauss_seidel = GaussSeidil(coefficients, constants, int(self.number_of_iterations.get()), eps,
                                           initial_array, float_converter)
                answer = gauss_seidel.solve()
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "Jacobi":
            try:
                jacobi = Jacobi(coefficients, constants, int(self.number_of_iterations.get()), eps, initial_array,
                                float_converter)
                answer = jacobi.solve()
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "Bisection":
            try:
                if self.bracketing_config.get() == "Xl & Xu":
                    bisection = Bisection(converter=float_converter,function=equation_str,rel_tolerance=eps, x_lower=float(self.x_lower.get()),
                                          x_upper=float(self.x_upper.get()), max_iterations=int(self.number_of_iterations.get()))
                the_plot = bisection.get_plot()
                answer = bisection.solve(0)
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "False-Position":
            try:
                if self.bracketing_config.get() == "Xl & Xu":
                    false_position = FalsePosition(converter=float_converter,function=equation_str,rel_tolerance=eps, x_lower=float(self.x_lower.get()),
                                          x_upper=float(self.x_upper.get()), max_iterations=int(self.number_of_iterations.get()))
                the_plot = false_position.get_plot()
                answer = false_position.solve(0)
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        elif self.method_selected.get() == "Fixed point":
            try:
                x = symbols('x')
                equation = simplify(equation_str)
                fixed_point = FixedPointIteration(equation, float(self.first_guess.get()), eps, int(self.number_of_iterations.get()), float_converter)
                answer, the_plot = fixed_point.fixedPointIteration()
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return
        elif self.method_selected.get() == "Newton-Raphson":
            try:
                x = symbols('x')
                equation = simplify(equation_str)
                newton_raphson = NewtonRaphson(equation, float(self.first_guess.get()), int(self.number_of_iterations.get()), eps, float_converter)
                answer, the_plot = newton_raphson.solve()
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return
        elif self.method_selected.get() == "Secant":
            try:
                x = symbols('x')
                equation = simplify(equation_str)
                secant = Secant(equation, float(self.first_guess.get()), float(self.second_guess.get()), int(self.number_of_iterations.get()), eps, float_converter)
                answer, the_plot = secant.solve()
            except ValueError:
                answer_window.destroy()
                error_window = tk.Toplevel(self)
                error_window.geometry("400x100")
                error_window.title("Error!")
                l1 = ttk.Label(error_window, text=f"Error!", font=14)
                l2 = ttk.Label(error_window, text=f"The method diverges!", font=12)
                l1.pack()
                l2.pack()
                b1 = ttk.Button(error_window, text="OK", command=close_error)
                b1.pack()
                return

        time_after = time.perf_counter_ns()
        time_taken = (time_after - time_before)

        l2 = ttk.Label(answer_window, text=f"Time Before = {time_before} ns")
        l2.grid(column=0, row=1)
        l3 = ttk.Label(answer_window, text=f"Time Taken = {time_taken} ns")
        l3.grid(column=3, row=1)
        l3 = ttk.Label(answer_window, text=f"Time After = {time_after} ns")
        l3.grid(column=6, row=1)
        results = ""
        if not (self.method_selected.get() in root_finding_methods):
            for i in range(len(answer)):
                results = results + f"{variables[i]} = {answer[i]}\n"
        else:
            plot_window = tk.Toplevel(self)
            plot_window.geometry("600x400")
            plot_window.title("The Plot")
            results = results + f"Root = {answer}"
            fig = the_plot
            canvas = FigureCanvasTkAgg(fig, master=plot_window)
            canvas.draw()
            canvas.get_tk_widget().pack()
            toolbar = NavigationToolbar2Tk(canvas, plot_window)
            toolbar.update()
            canvas.get_tk_widget().pack()

        l4 = ttk.Label(answer_window, text=results, font=10)
        l4.grid(column=3, row=2)

        b1 = ttk.Button(answer_window, text="OK", command=close_answer)
        b1.grid(column=3, row=5)


if __name__ == "__main__":
    app = App()
    app.mainloop()

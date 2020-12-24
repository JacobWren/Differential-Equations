# Plotting Direction Fields to investigate the behavior of solutions without solving autonomous differential equations.
# Actually, we do solve the ODE's since we plot a few projectories (actual solutions).

import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction
from scipy.integrate import odeint


# Obtain user input.
correct_de = input("Is your differential equation of the form: dv/dt = a + b*v for some real numbers a & b? "
                         "(answer: 'Yes' or 'No' ): ")
if correct_de == "No":
    print("Sorry, this program will not work for your differential equation.")
elif correct_de == "Yes":
    a = input("Enter a: ")
    if '/' in a:
        a = float(Fraction(a))
    else:
        a = float(a)
    print()
    b = input("Enter b: ")
    if '/' in b:
        b = float(Fraction(b))
    else:
        b = float(b)
    print()

# Retry... if necessary.
while correct_de != "No" and correct_de != "no" and correct_de != "Yes" and correct_de != "yes":
    correct_de = input("Please enter your response again: ")


def Equilibrium():
    # Find equilibrium solution (i.e., terminal velocity),
    # in other words what value of v will cause dv/dt to be zero?
    Eq_Sol = a / -b
    return Eq_Sol


# The equilibrium solution will dictate the plot range.
E = Equilibrium()
y_lower = E * .70
y_upper = E * 1.30

x_lower = 0
x_upper = max(int((1 / (y_upper - y_lower)) * 10), 10) # If the range is large, the domain should be smaller, but never
# less than 10.


def Tangent_finder():
    y_prime = [] # Initialize list that will contain the slope of a solution for a given y.
    line_info = [] # We will store information of the tangent line here.

    for y in np.linspace(y_lower, y_upper, 10): # This ensures we will have
        # 20 points in any vertical "column".
        slope = a + (b * y) # Find slope.
        y_prime.append(slope)
        # x range (x and t are used synonymously)
        for t in np.linspace(x_lower, x_upper, (x_upper - x_lower)): # This ensures we will have '(x_upper - x_lower)' "columns".
            y_int = y - slope * t # Need to find y-intercept for each "line" (slope) that we will plot.
            line_info.append([y, slope, t, y_int]) # Contains y-value, slope, x-value, and y-intercept.
    return line_info


def Plot_tangents():
    line_info = Tangent_finder() # Make call to Tangent_finder function.

    # We need to plot many "mini" lines, i.e., slopes.
    # Remember these lines are tangent to the solution at 'y' (the first argument of line_info),
    # so the lines we plot should be centered around that 'y' value

    for line in line_info: # We will plot every line computed above...
        # to do this we will find 2 points for each line.
        # This determines the length of the vectors.
        x_distance = min(.3, abs(.08 * (y_upper - y_lower))) # We need this to be positive in order to place the arrow on the
        # correct end. In other words, we need x_cord2 > x_cord1.
        # Here, we are trying to center the line around the 'y' value.
        x_cord1 = line[2] - x_distance
        x_cord2 = line[2] + x_distance

        # This is just for clarity.
        slope = line[1]
        y_int = line[3]

        # Compute y values
        y_cord1 = x_cord1 * slope + y_int # y = mx + b!
        y_cord2 = x_cord2 * slope + y_int

        # Plot our 2 points computed above.
        plt.plot([x_cord1, x_cord2], [y_cord1, y_cord2])

        # Let's add an arrow to the tip of the slope.
        # For arrow placement we need to specify a delta x and a delta y.
        delta_x = .001
        delta_y = y_cord2 - ((x_cord2 - delta_x) * slope + y_int) # Change in y-direction.
        # 'xytext' is the starting direction of the arrow.
        # A negative 'headlength' refelcts the arrow inward.
        plt.annotate("", xy = (x_cord2 + delta_x, y_cord2 + delta_y), xytext = (x_cord2, y_cord2),
                     arrowprops = dict(headwidth = 3, headlength = min(abs(y_upper - y_lower) * 2, 6), width = (y_upper - y_lower) * 2))

    plt.xlabel("t",
             fontweight = 'bold',
             fontsize = 14,
               loc = "right") # 't' for time.
    plt.ylabel("v",
             fontweight = 'bold',
             rotation = 0,
             fontsize = 14,
             loc= "top") # 'v' for velocity.

    plt.xlim(x_lower, x_upper)
    plt.ylim(y_lower, y_upper)
    plt.margins(x = 0) # Removes horizontal spacing from plot.

    # Plot trajectories (solutions).
    # Initial conditions.
    if b < 0:
        v0_0 = np.linspace(1.05 * E, .95 * E, 2, endpoint=True)
        v0_1 = np.linspace(1.1 * E, .9 * E, 2, endpoint=True)
        v0_2 = np.linspace(1.15 * E, .85 * E, 2, endpoint=True)
        v0_3 = np.linspace(y_lower, y_upper, 3, endpoint=True)
        v0_4 = np.linspace(y_lower * -2, y_upper * 2, 3, endpoint=True)
        v0_5 = np.linspace(y_lower * -5, y_upper * 5, 3, endpoint=True)
        v0_6 = np.linspace(y_lower * -10, y_upper * 10, 3, endpoint=True)
        v0_7 = np.linspace(y_lower * -100, y_upper * 100, 3, endpoint=True)
        v0_8 = np.linspace(y_lower * -1000, y_upper * 1000, 3, endpoint=True)
        v0_9 = np.linspace(y_lower * -10000, y_upper * 10000, 3, endpoint=True)
        v0_10 = np.linspace(y_lower * -100000, y_upper * 100000, 3, endpoint=True)
        v0_11 = np.linspace(y_lower * -1000000, y_upper * 1000000, 3, endpoint=True)
        v0_12 = np.linspace(y_lower * -10000000, y_upper * 10000000, 3, endpoint=True)

        I_C_list = [v0_0, v0_1, v0_2, v0_3, v0_4, v0_5, v0_6, v0_7, v0_8, v0_9, v0_10, v0_11, v0_12]

        # We need this outer for loop to space out the initial conditions, there may be a better way to do this...
        # Each run through the outer loop, uses a different array of initial conditions, where the solutions tend to be
        # clustered together for each array - hence why we need several different arrays with different ranges for
        # proper spacing.
        for I_Cs in I_C_list:
            for I_C in I_Cs:  # We didn't specify a particular initial condition,
                # rather we will look at solutions to many such I.C.'s
                # This function is required for odeint().
                def ODE(v, t):
                    dvdt = a + b * v
                    return dvdt

                t = np.linspace(x_lower, x_upper)
                # Find v for a particular I.C.
                v = odeint(ODE, I_C, t)
                plt.plot(t, v)
    else: # b positive case
        l_1 = [1.0000000001, 1.000000001, 1.00000001, 1.0000001, 1.000001, 1.00001, 1.0001, 1.001, 1.01]
        v0_0 = np.linspace(l_1[0] * E, (1 - (l_1[0] - 1)) * E, 2, endpoint=True)
        v0_1 = np.linspace(l_1[1] * E, (1 - (l_1[1] - 1)) * E, 2, endpoint=True)
        v0_2 = np.linspace(l_1[2] * E, (1 - (l_1[2] - 1)) * E, 2, endpoint=True)
        v0_3 = np.linspace(l_1[3] * E, (1 - (l_1[3] - 1)) * E, 2, endpoint=True)
        v0_4 = np.linspace(l_1[4] * E, (1 - (l_1[4] - 1)) * E, 2, endpoint=True)
        v0_5 = np.linspace(l_1[5] * E, (1 - (l_1[5] - 1)) * E, 2, endpoint=True)
        v0_6 = np.linspace(l_1[6] * E, (1 - (l_1[6] - 1)) * E, 2, endpoint=True)
        v0_7 = np.linspace(l_1[7] * E, (1 - (l_1[7] - 1)) * E, 2, endpoint=True)
        v0_8 = np.linspace(l_1[8] * E, (1 - (l_1[8] - 1)) * E, 2, endpoint=True)

        I_C_list = [v0_0, v0_1, v0_2, v0_3, v0_4, v0_5, v0_6, v0_7, v0_8]

        for I_Cs in I_C_list:
            for I_C in I_Cs:
                def ODE(v, t): # Even though the ODE is autonomous, t is needed.
                    dvdt = a + b * v
                    return dvdt

                t = np.linspace(x_lower, x_upper)
                v = odeint(ODE, I_C, t)
                plt.plot(t, v)

    locs, labels = plt.yticks()  # locs gives us the y-tick values on the plot.
    plt.axhline(E, x_lower, x_upper, linewidth = 2, color = 'black', label = 'Equilibrium Solution') # Plot equilibrium solution.
    label = True
    for tick in locs:
        if tick < 1.03 * E and tick > .97 * E:
            label = False
    if label: # We don't need to label the equilibrium value if the plot will have a tick mark close by.
        plt.text(-.71, E, str(E)) # Here we append the equilibrium value to the plot.
    plt.legend(loc = 'lower right', bbox_to_anchor = (.75, -.155)) 
    plt.title('Direction field', fontweight = 'bold',
              fontsize = 18)
    plt.show()

Plot_tangents()
import math as m
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from numpy import random
from orbital_functions import orbit, get_shortest, n_orbit_positions, get_coordinates, read_input_list, get_moid, fig, ax
import pprint
import sys
import time

plt.title("Calculte Earth MOID using sliders . . .")

# initial variables
e = .01
a = 0 
loan = 0
aop = 0
i = 0
opacity = 0.75
color = 'black'

# create slider widgets
slider_area = fig.subplots_adjust(left = 0.1, bottom = 0.3)    
ax_slider1 = plt.axes([0.05, 0.55, 0.2, 0.025]) # arguments: x, y, width, height
ax_slider2 = plt.axes([0.05, 0.5, 0.2, 0.025]) # og: [0.1, 0.1, 0.8, 0.05] -- > (thinner) [0.325, 0.1, 0.35, 0.05]
ax_slider3 = plt.axes([0.05, 0.45, 0.2, 0.025])

ax_slider4 = plt.axes([0.05, 0.65, 0.2, 0.025])
ax_slider5 = plt.axes([0.05, 0.6, 0.2, 0.025])

# look into using text field input for these variables instead . . .
i_slider = Slider(ax_slider1, 'i', valmin = 0, valmax = 180, valinit = 0, valstep = 1) # 10
loan_slider = Slider(ax_slider2, 'loan', valmin = 0, valmax = 360, valinit = 0, valstep = 1) # 20
aop_slider = Slider(ax_slider3, 'aop', valmin = 0, valmax = 360, valinit = 0, valstep = 1) # 20
a_slider = Slider(ax_slider4, 'a', valmin = 0, valmax = 3, valinit = 0, valstep = 0.01) # 0.1 --> max : 10
e_slider = Slider(ax_slider5, 'e', valmin = 0.01, valmax = 0.99, valinit = 0.01, valstep = 0.001) # 0.01

# earth
earth_orbit = orbit('Earth', 102.94719, -11.26064, 0.00005, 0.01671022, 1.00000011, n_orbit_positions, 'blue', True, False, 1.0)

# hypothetical orbit
hypothetical_orbit = orbit('orbit', aop, loan, i, e, a, n_orbit_positions, color, True, False, opacity)

def update_axis():
    Axes3D.cla(ax)

    ax.set_xlim3d(-2.5, 2.5)
    ax.set_ylim3d(-2.5, 2.5)
    ax.set_zlim3d(-2.5, 2.5)

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

# orbit(name, loan, aop, i, e, a, n, color, plot?, print?, opacity)
def i_update(n):
    a, e, loan, aop = a_slider.val, e_slider.val, loan_slider.val, aop_slider.val

    update_axis()

    updated_earth_orbit = orbit('Earth', 102.94719, -11.26064, 0.00005, 0.01671022, 1.00000011, n_orbit_positions, 'blue', True, False, 1.0)
    updated_hypothetical_orbit = orbit('orbit', aop, loan, n, e, a, n_orbit_positions, color, True, False, 0.75)
    updated_moid = get_moid(earth_orbit, updated_hypothetical_orbit, 'Earth', 'hypothetical_orbit', True, False)

    ax.set_title("Earth moid of orbit: " + str(updated_moid))

def loan_update(n):
    a, e, i, aop = a_slider.val, e_slider.val, i_slider.val, aop_slider.val
    
    update_axis()

    updated_earth_orbit = orbit('Earth', 102.94719, -11.26064, 0.00005, 0.01671022, 1.00000011, n_orbit_positions, 'blue', True, False, 1.0)
    updated_hypothetical_orbit = orbit('orbit', aop, n, i, e, a, n_orbit_positions, color, True, False, 0.75)
    updated_moid = get_moid(earth_orbit, updated_hypothetical_orbit, 'Earth', 'hypothetical_orbit', True, False)

    ax.set_title("Earth moid of orbit: " + str(updated_moid))

def aop_update(n):
    a, e, i, loan = a_slider.val, e_slider.val, i_slider.val, loan_slider.val
    
    update_axis()

    updated_earth_orbit = orbit('Earth', 102.94719, -11.26064, 0.00005, 0.01671022, 1.00000011, n_orbit_positions, 'blue', True, False, 1.0)
    updated_hypothetical_orbit = orbit('orbit', n, loan, i, e, a, n_orbit_positions, color, True, False, 0.75)
    updated_moid = get_moid(earth_orbit, updated_hypothetical_orbit, 'Earth', 'hypothetical_orbit', True, False)

    ax.set_title("Earth moid of orbit: " + str(updated_moid))

def a_update(n):
    e, i, loan, aop = e_slider.val, i_slider.val, loan_slider.val, aop_slider.val

    update_axis()

    updated_earth_orbit = orbit('Earth', 102.94719, -11.26064, 0.00005, 0.01671022, 1.00000011, n_orbit_positions, 'blue', True, False, 1.0)
    updated_hypothetical_orbit = orbit('orbit', aop, loan, i, e, n, n_orbit_positions, color, True, False, 0.75)
    updated_moid = get_moid(earth_orbit, updated_hypothetical_orbit, 'Earth', 'hypothetical_orbit', True, False)

    ax.set_title("Earth moid of orbit: " + str(updated_moid))

def e_update(n):
    a, i, loan, aop = a_slider.val, i_slider.val, loan_slider.val, aop_slider.val
    
    update_axis()

    updated_earth_orbit = orbit('Earth', 102.94719, -11.26064, 0.00005, 0.01671022, 1.00000011, n_orbit_positions, 'blue', True, False, 1.0)
    updated_hypothetical_orbit = orbit('orbit', aop, loan, i, n, a, n_orbit_positions, color, True, False, 0.75)
    updated_moid = get_moid(earth_orbit, updated_hypothetical_orbit, 'Earth', 'hypothetical_orbit', True, False)

    ax.set_title("Earth moid of orbit: " + str(updated_moid))

i_slider.on_changed(i_update)
loan_slider.on_changed(loan_update)
aop_slider.on_changed(aop_update)
a_slider.on_changed(a_update)
e_slider.on_changed(e_update)

plt.show()
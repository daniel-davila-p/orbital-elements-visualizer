'''
 Program that calculates an object's MOID with Jupiter.
 MOID = minimum orbit intersection distance. ie, the shortest 
 distance between the object's orbit and Jupiter's orbit.

 inputs: a    - semimajor axis in AU
         e    - eccentricity
         i    - inclination in degrees
         aop  - argument of perihelion in degrees
         loan - longitude of ascending node in degrees
         jloc - 3 x N matrix of (xyz) positions of Jupiter's 
                orbit at N different locations

 outputs: oid - matrix that holds all pairwise distances between
          a point on the object's orbit and a point on Jupiter's 
          orbit.

          MOID - the minimum of the oid matrix

 The conversion from orbital elements to Cartesian positions is done
 like this: It the object's own frame, consider that its orbit is in
 its own xy plane. Use true anomaly f and heliocentric distance r
 as the polar coordiantes. True anomaly f runs from 0 to 360 deg.

 From f you can calculate eccentric anomaly E, and from that you
 can get helioc. distance r.

 Note that f is essentially the angular coordinate in the plane of
 the object's orbit. In other words, if we were to describe the position
 of a location in the object's orbit in spherical coordinates, we could
 use r = heliocentric distance as the radial coordinate, and f as the
 angular coordinates in the orbit plane.
'''
import sys
import math as m 
import numpy as np 
import pprint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime

n_orbit_positions = 1000
pi = m.pi

fig = plt.figure(figsize=(6, 3))
ax = plt.axes(projection="3d")
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

ax.set_xlim3d(-2.5, 2.5)
ax.set_ylim3d(-2.5, 2.5)
ax.set_zlim3d(-2.5, 2.5)

ax.w_xaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
ax.w_yaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))
ax.w_zaxis.set_pane_color((0.8, 0.8, 0.8, 1.0))

def plot_body(x, y, z, color):
    u = np.linspace(0, 2 * pi, 100)
    v = np.linspace(0, pi, 100)

    body_x = 0.05 * np.outer(np.cos(u), np.sin(v))
    body_y = 0.05  * np.outer(np.sin(u), np.sin(v))
    body_z = 0.025 * np.outer(np.ones(np.size(u)), np.cos(v))
        
    ax.plot_surface(body_x + x, body_y + y, body_z + z, rstride=4, cstride=4, color=color)

def bubble_sort(array):
    size = len(array)

    for i in range(size - 1):
        for j in range(0, size - i - 1):
            if array[j] > array[j + 1]:
                array[j], array[j + 1] = array[j + 1], array[j]

def get_distance(x1, y1, z1, x2, y2, z2):
    distance = m.sqrt(m.pow(x2 - x1, 2) + m.pow(y2 - y1, 2) + m.pow(z2 - z1, 2))
    return distance

def get_shortest(list1, list2, name1, name2, opacity, plot_variable, print_text):
    shortest = sys.maxsize
 
    for i in range(0, len(list1)):
        for j in range(0, len(list2)):
 
            x1 = list1[i][0]
            y1 = list1[i][1]
            z1 = list1[i][2]
 
            x2 = list2[j][0]
            y2 = list2[j][1]
            z2 = list2[j][2]

            current_distance = get_distance(x1, y1, z1, x2, y2, z2)

            if current_distance < shortest:
                shortest = current_distance
                x_line = x1, x2
                y_line = y1, y2
                z_line = z1, z2

    if plot_variable == True:
        ax.plot3D(x_line, y_line, z_line, 'red', alpha = opacity)
    if print_text == True:
        print("After " + str(m.pow(n_orbit_positions, 2)) + " comparisons . . . " + name1 + " MOID " + name2 + " " + str(round(shortest, 3)) + " when:")
        print("\t\t\t\t\t" + name1 + ": (" + str(x1) + ", " + str(y1) + ", " + str(z1) + ")")
        print("\t\t\t\t\t" + name2 + ": (" + str(x2) + ", " + str(y2) + ", " + str(z2) + ")")
    
    return shortest

# returns set of coordinates (0, 2pi) option to plot, decide it's transparency, and color
def orbit(name, loan, aop, i, e, a, n, color, plot_variable, print_to_file, opacity):
    if print_to_file:
        ifp = open(name + "_orbit_NO_DATE.txt", "w")
    
    # make some sines and cosines of the angle orbital elements
    cloan = m.cos((loan * pi) / 180)
    sloan = m.sin((loan * pi) / 180)

    caop = m.cos((aop * pi) / 180)
    saop = m.sin((aop * pi) / 180)

    ci = m.cos((i * pi) / 180)
    si = m.sin((i * pi) / 180)

    # now make the rotation matrix
    '''
    # Murray and Dermott book   
    coeff_x = [(cloan * caop) - (ci * sloan * saop), (cloan * saop) - (ci * sloan * caop), (si * sloan)]
    coeff_y = [(sloan * caop) + (ci * cloan * saop), (sloan * saop) + (ci * cloan * caop), (si * cloan)]
    coeff_z = [(si * saop), (si * caop), ci]
    '''
    # wikipedia version
    coeff_x = (cloan * caop) - (ci * sloan * saop), (sloan * caop) + (ci * cloan * saop), (si * saop)
    coeff_y = (-cloan * saop) - (ci * sloan * caop), (-sloan * saop) + (ci * cloan * caop), (si * caop)
    coeff_z = (si * sloan), (-si * cloan), ci
    
    f = []
    sub1 = []
    sub2 = []
    tan_e_2 = []
    ecc_anom = []
    r = []

    for i in range(0, n): # f is the True Anomaly
        current_f = (i / n) * 360
        current_sub1 = m.tan(current_f * (pi / 180) / 2)
        current_sub2 = m.sqrt((1 + e) / (1 - e))
        current_tan_e_2 = current_sub1 / current_sub2
        current_ecc_anom = 2 * m.atan(current_tan_e_2)
        current_r = a * (1 - (e * m.cos(current_ecc_anom)))

        f.append(current_f)
        sub1.append(current_sub1)
        sub2.append(current_sub2)
        tan_e_2.append(current_tan_e_2)
        ecc_anom.append(current_ecc_anom)
        r.append(current_r)# spherical coordinates created

    pos_array = np.empty(shape=(3, n_orbit_positions),dtype='object') 

    for i in range(0, n_orbit_positions):
        pos_array[0][i] = (r[i] * m.cos(f[i] * (pi / 180)))
        pos_array[1][i] = (r[i] * m.sin(f[i] * (pi / 180))) 
        pos_array[2][i] = 0

    x = np.matmul(coeff_x, pos_array)
    y = np.matmul(coeff_y, pos_array)
    z = np.matmul(coeff_z, pos_array)

    coordinates = []
  
    for i in range(0, n_orbit_positions):
        current_coordinate = round(x[i], 5), round(y[i], 5), round(z[i], 5) #str(date)]
        coordinates.append(current_coordinate)

    if plot_variable:
        if print_to_file:
            ifp.write(name + " orbit\n********************************************\n\n")
            for i in range(0, n_orbit_positions):
                ifp.write(str(coordinates[i]) + "\n") 
        ax.plot3D(x, y, z, color, alpha = opacity)

        #Axes3D.cla(ax)
        
        #plt.show()
    if print_to_file:
        ifp.close()

    return coordinates

'''
 Now that we have a set of coordinates, we need to determine what time these orbits start, from JPL Horizons 
 parameters: 
            - seconds from 2021-01-01 00:00:00
            - end time from 2021-01-01 00:00:00
            - orbital speed in au/s
'''
def get_coordinates(name, list, orbital_speed, start_time, end_time):
    ifp = open(name + "_orbit_WITH_DATE", "w")
    total_distance = 0
    i = 0 
    initial = 0
    d_time = 0
    new_list = []

    while start_time + d_time < end_time:
        if i == len(list):
            i = 0
            initial = 0
        if i > 1:
            initial += 1   

        total_distance = get_distance(list[initial][0], list[initial][1], list[initial][2], list[i][0], list[i][1], list[i][2]) # in au
        #print(total_distance)
        d_time += (total_distance / orbital_speed) # in seconds
        date = datetime.fromtimestamp(start_time + d_time).date() # time since 2021-01-01 00:00:00

        new_coordinate = [list[initial][0], list[initial][1], list[initial][2], str(date)]
        new_list.append(new_coordinate)
        i += 1

    for i in range(0, len(new_list)):
        ifp.write(str(new_list[i]) + "\n") 

    return new_list

# read in the positional data from the input file
def read_input_list(input_file_name):
    input_file = open(input_file_name, "r")
    list_coords = []
    previous_d = 0
    for line in input_file:
        row = line.split(",")
        current_x, current_y, current_z, current_d = row[0].replace("[", ""), float(row[1]), float(row[2]), row[3].replace("]", "")
        if current_d != previous_d:
            current_coordinate = float(current_x), float(current_y), float(current_z), current_d
            list_coords.append(current_coordinate)
        previous_d = current_d
    #pprint.pprint(list_coords)
    return list_coords

def read_input_list_contour(input_file_name):
    input_file = open(input_file_name, "r")
    list_coords = []
    
    for line in input_file:
        row = line.split(",")
        current_moid, current_i, current_loan, current_aop = float(row[0].replace("(", "")), float(row[1]), float(row[2]), float(row[3].replace(")", ""))
        
        current_coordinate = current_moid, current_i, current_loan, current_aop
        list_coords.append(current_coordinate)

    plt.close()
    return list_coords

def read_input_list_new_contour(input_file_name):
    input_file = open(input_file_name, "r")
    list_coords = []
    
    for line in input_file:
        row = line.split(",")
        current_moid, current_i, current_aop, current_a, current_e = float(row[0].replace("(", "")), float(row[1]), float(row[2]), float(row[3]), float(row[4].replace(")", ""))
        
        current_coordinate = current_moid, current_i, current_aop, current_a, current_e
        list_coords.append(current_coordinate)

    plt.close()
    return list_coords

def get_moid(list1, list2, name1, name2, plot_variable, print_text):
    moid = sys.maxsize
 
    for i in range(0, len(list1)):
        for j in range(0, len(list2)):
 
            x1, y1, z1 = list1[i][0], list1[i][1], list1[i][2]
            x2, y2, z2 = list2[j][0], list2[j][1], list2[j][2]
            
            current_distance = get_distance(x1, y1, z1, x2, y2, z2)

            if current_distance < moid:
                moid = current_distance
                x_line = x1, x2
                y_line = y1, y2
                z_line = z1, z2

                body1_moidx, body1_moidy, body1_moidz = x1, y1, z1
                body2_moidx, body2_moidy, body2_moidz = x2, y2, z2

    if plot_variable == True:
        ax.plot3D(x_line, y_line, z_line, 'red', alpha = 1.0)
        # plot_body(body1_moidx, body1_moidy, body1_moidz, 'black')
        # plot_body(body2_moidx, body2_moidy, body2_moidz, 'black')
        '''
        ax.text(body1_moidx * 1.5, body1_moidy, body1_moidz, name1, color = 'black', size = 10)
        ax.text(body2_moidx * 1.5, body2_moidy, body2_moidz, name2, color = 'black', size = 10)
        ax.text2D(1, 1, "OID of " + str(moid), transform=ax.transAxes)
        ax.text2D(1, 0.85, "reached on " + str(shortest_date), transform=ax.transAxes)
        '''

    if print_text == True:
        print("MOID " + name1 + " & " + name2 + " between " + str(list1[0][3].strip("\n")) + " and " + str(list1[len(list1) - 1][3]).strip("\n") + ": " + str(moid))
        # print("Reached on " + shortest_date)
    return moid


# function for plotting vector
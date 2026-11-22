import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import quad

# (typical) Coefficients of drag for varius shapes
# One of these coefficients will be used for each profile based on the shape it most clearly resembles
cd_sphere = .47
cd_cube = 1.05
cd_cylinder = .82
cd_cone = .50
cd_plate = 1.28
cd_half_sphere = .42

def calculate_drag_force(v_x, v_y, rho, A_front, A_bottom, Cdf, Cdb):
    #Calculates the speed vector magnitude
    v_mag = math.hypot(v_x, v_y)

    #Avoids division by zero if the objecft is no longer moving
    if v_mag == 0:
        return 0.0, 0.0

    #Calculates the direction cosines to scale the front and bottom areas
    cos = abs(v_x) / v_mag
    sin = abs(v_y) / v_mag

    #Calculates the effective drag area
    CdA = (Cdf * A_front * cos) + (Cdb * A_bottom * sin)

    #Calculates the magnitude of the total drag force
    F_drag = .5 * rho * CdA * (v_mag ** 2)

    #Projects the total drag force onto the x and y axis and reverses the direction
    F_drag_x = -F_drag * (v_x / v_mag)
    F_drag_y = -F_drag * (v_y / v_mag)

    return F_drag_x, F_drag_y

def generate_random_profile(length, num_points=6):
    """
    Generates a random aerodynamic profile of a given length
    """

    # Creates evenly spaced x coordinates from 0 to length
    x_coords = np.linspace(0, length, num_points)

    # Generates random y coordinates for the profile (thickness)
    y_coords = np.random.uniform(0.2, 2.5, num_points)

    # Makes the front and back of the profile to be at y=0
    y_coords[0] = 0
    y_coords[-1] = 0

    # Creates a continuous function through the points to represent the shape with smooth ends
    random_shape = CubicSpline(x_coords, y_coords, bc_type='clamped')

    # How to use / Example: random_side_profile = generate_random_profile(length=10)
    # Random profile that is 10 meters long

    return random_shape


# Potential function for calculating the drag coefficient

def estimate_cd(profile, length, area):
    """
    Estimates the drag coefficient of a given profile based on its shape and area
    """

    #Take 200 evenly spaced points as samples along length and finds peak half-width
    x_values = np.linspace(0, length, 200)
    max_y = max(abs(profile(x_values)))
    max_thickness = 2 * max_y #Profile is symmetric about the centerline

    #Gives minimum Cd if the profile is completely flat
    if max_thickness == 0:
        return .05

    #Computes bounding rectangle area enclosing full profile
    bound_area = length * max_thickness
    fullness = area / bound_area #Measures the bluntness / bloxiness of the shpae
    thickness_ratio = max_thickness / length

    cd_est = .15 + (.75 * fullness) + (.2 * min(thickness_ratio, 1))
    
    return cd_est


# Example Run putting everything together

length = 2.0

front_profile = generate_random_profile(length)
bottom_profile = generate_random_profile(length)

A_front = 2 * quad(front_profile, 0, length)[0]  # Area under the front profile curve
A_bottom = 2 * quad(bottom_profile, 0, length)[0]  # Area under the bottom profile curve

Cdf = estimate_cd(front_profile, length, A_front)
Cdb = estimate_cd(bottom_profile, length, A_bottom)

dx, dy = calculate_drag_force(v_x=10, v_y=5, rho=1.225, A_front=A_front, A_bottom=A_bottom, Cdf=Cdf, Cdb=Cdb)

print(f"Drag force in x-direction: {dx:.2f} N")
print(f"Drag force in y-direction: {dy:.2f} N")

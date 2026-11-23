import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import quad

# Constants
rho = 1.225
g = 9.807


def calculate_drag_force(v_x, v_y, rho, A_front, A_bottom, Cdf, Cdb):
    # Calculates the speed vector magnitude
    v_mag = math.hypot(v_x, v_y)

    # Avoids division by zero if the object is no longer moving
    if v_mag == 0:
        return 0.0, 0.0

    # Calculates the direction cosines to scale the front and bottom areas
    cos = abs(v_x) / v_mag
    sin = abs(v_y) / v_mag

    # Calculates the effective drag area
    CdA = (Cdf * A_front * cos) + (Cdb * A_bottom * sin)

    # Calculates the magnitude of the total drag force
    F_drag = 0.5 * rho * CdA * (v_mag ** 2)

    # Projects the total drag force onto the x and y axis and reverses the direction
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

    return random_shape


def estimate_cd(profile, length, area):
    """
    Estimates the drag coefficient of a given profile based on its shape and area
    """
    # Take 200 evenly spaced points as samples along length and finds peak half-width
    x_values = np.linspace(0, length, 200)
    max_y = max(abs(profile(x_values)))
    max_thickness = 2 * max_y  # Profile is symmetric about the centerline

    # Gives minimum Cd if the profile is completely flat
    if max_thickness == 0:
        return 0.05

    # Computes bounding rectangle area enclosing full profile
    bound_area = length * max_thickness
    fullness = area / bound_area  # Measures the bluntness / boxiness of the shape
    thickness_ratio = max_thickness / length

    cd_est = 0.15 + (0.75 * fullness) + (0.2 * min(thickness_ratio, 1))
    return cd_est


# EULER SIMULATION FUNCTION
def trajectory_euler(p0=(0, 0), v0=(10, 10), mass=3, A_front=1.0, A_bottom=1.0,
                     Cdf=0.47, Cdb=0.47, h=0.01, stop_height=0, max_steps=100000, plot=True):
    # Validate inputs
    if h <= 0:
        raise ValueError("h must be positive")
    if mass <= 0:
        raise ValueError("Mass must be positive")

    # Initial state [x, y, vx, vy]
    state = np.array([*p0, *v0], dtype=float)
    positions = [state[:2].copy()]
    velocities = [state[2:].copy()]
    times = [0.0]

    # Derivative Function
    def deriv(state):
        x, y, vx, vy = state
        F_drag_x, F_drag_y = calculate_drag_force(vx, vy, rho, A_front, A_bottom, Cdf, Cdb)

        ax = F_drag_x / mass
        ay = -g + (F_drag_y / mass)

        return np.array([vx, vy, ax, ay])

    # EULER'S METHOD
    t = 0.0
    for i in range(max_steps):
        previous = state.copy()

        # Euler integration
        state = state + h * deriv(state)
        t_next = t + h

        # Ground collision
        if state[1] <= stop_height:
            denom = previous[1] - state[1]
            alpha = (
                (previous[1] - stop_height) / denom if denom != 0 else 0
            )

            hit_state = previous + alpha * (state - previous)
            hit_state[1] = stop_height

            positions.append(hit_state[:2].copy())
            velocities.append(hit_state[2:].copy())
            times.append(t + alpha * h)
            break

        positions.append(state[:2].copy())
        velocities.append(state[2:].copy())
        times.append(t_next)
        t = t_next
    else:
        raise RuntimeError("Maximum steps exceeded")

    # CONVERT RESULTS INTO TABLE
    p = np.array(positions)
    v = np.array(velocities)
    t = np.array(times)

    table = pd.DataFrame({
        "Time (s)": t,
        "X Position (m)": p[:, 0],
        "Y Position (m)": p[:, 1],
        "X Velocity (m/s)": v[:, 0],
        "Y Velocity (m/s)": v[:, 1]
    })

    # GRAPH RESULTS
    if plot:
        # Graph 1: Trajectory
        plt.figure(figsize=(10, 6))
        plt.plot(p[:, 0], p[:, 1])
        plt.title("Projectile Trajectory - Euler's Method")
        plt.xlabel("Horizontal Distance (m)")
        plt.ylabel("Height (m)")
        plt.grid()
        plt.show()

        # Graph 2: Velocity
        plt.figure(figsize=(10, 6))
        plt.plot(t, v[:, 0], label="Horizontal Velocity")
        plt.plot(t, v[:, 1], label="Vertical Velocity")
        plt.title("Velocity vs Time - Euler's Method")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity (m/s)")
        plt.legend()
        plt.grid()
        plt.show()

    # RETURN RESULTS TABLE
    return table


# Generate Random Profile & Aerodynamic Parameters
length = np.random.uniform(0.5, 3.0)
front_profile = generate_random_profile(length)
bottom_profile = generate_random_profile(length)

A_front = 2 * quad(front_profile, 0, length)[0]
A_bottom = 2 * quad(bottom_profile, 0, length)[0]

Cdf = estimate_cd(front_profile, length, A_front)
Cdb = estimate_cd(bottom_profile, length, A_bottom)

# RUN EULER SIMULATION
euler_trajectory = trajectory_euler(
    p0=(0, 0),
    v0=(10, 10),
    mass=3,
    A_front=A_front,
    A_bottom=A_bottom,
    Cdf=Cdf,
    Cdb=Cdb,
    h=0.01
)
print(euler_trajectory)

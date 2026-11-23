
import numpy as np
import matplotlib.pyplot as plt

# Constants

rho = 1.225
g = 9.807
Cd = 0.47


# Euler's Method

def Euler(state, h, deriv):

    return state + h * deriv(state)


# Projectile Simulation

def trajectory_sphere(p0, v0, mass, radius, h,
                      stop_height=0, max_steps=100000):

    if h <= 0:
        raise ValueError("h must be positive")
    if mass <= 0:
        raise ValueError("Mass must be positive")
    if radius < 0:
        raise ValueError("Radius cannot be negative")

    area = np.pi * radius**2

    # Initial state: [x, y, vx, vy]
    state = np.array([*p0, *v0], dtype=float)

    positions = [state[:2].copy()]
    velocities = [state[2:].copy()]
    times = [0.0]

    # Derivative Function

    def deriv(state):

        x, y, vx, vy = state

        speed = np.hypot(vx, vy)

        # Drag acceleration factor
        drag_factor = 0.5 * Cd * rho * area * speed / mass

        ax = -drag_factor * vx
        ay = -g - drag_factor * vy

        # [dx/dt, dy/dt, dvx/dt, dvy/dt]
        return np.array([vx, vy, ax, ay])

    # Integration Loop

    t = 0.0

    for i in range(max_steps):

        previous = state.copy()

        # Advance simulation using Euler's Method
        state = Euler(state, h, deriv)

        t_next = t + h

        # Check ground collision
        if state[1] <= stop_height:

            denom = previous[1] - state[1]

            alpha = (
                (previous[1] - stop_height) / denom
                if denom != 0 else 0
            )

            # Interpolate state at impact
            hit_state = previous + alpha * (state - previous)
            hit_state[1] = stop_height

            positions.append(hit_state[:2].copy())
            velocities.append(hit_state[2:].copy())
            times.append(t + alpha * h)

            return positions, velocities, times

        # Record current state
        positions.append(state[:2].copy())
        velocities.append(state[2:].copy())
        times.append(t_next)

        t = t_next

    raise RuntimeError("Maximum steps exceeded")


# Run Simulation

p, v, t = trajectory_sphere(
    p0=(0, 0),
    v0=(10, 10),
    mass=3,
    radius=0.5,
    h=0.01
)


# Get results

p = np.array(p)
v = np.array(v)
t = np.array(t)


# Graph 1: Trajectory

plt.figure(figsize=(10, 6))

plt.plot(p[:, 0], p[:, 1])

plt.title("Sphere Trajectory - Euler's Method")
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



# new code to help make easier for implementing multiple shapes

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# SPHERE SIMULATION FUNCTION
def sphere(p0=(0, 0), v0=(10, 10), mass=3, radius=0.5,
           h=0.01, stop_height=0, max_steps=100000,
           plot=True):

    # Constants
    rho = 1.225
    g = 9.807
    Cd = 0.47

    # Validate inputs
    if h <= 0:
        raise ValueError("h must be positive")
    if mass <= 0:
        raise ValueError("Mass must be positive")
    if radius < 0:
        raise ValueError("Radius cannot be negative")

    # Cross-sectional area of sphere
    area = np.pi * radius**2

    # Initial state [x, y, vx, vy]
    state = np.array([*p0, *v0], dtype=float)

    positions = [state[:2].copy()]
    velocities = [state[2:].copy()]
    times = [0.0]

    # Derivative Function
    def deriv(state):

        x, y, vx, vy = state

        speed = np.hypot(vx, vy)

        # Drag acceleration factor
        drag_factor = 0.5 * Cd * rho * area * speed / mass

        ax = -drag_factor * vx
        ay = -g - drag_factor * vy

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
                (previous[1] - stop_height) / denom
                if denom != 0 else 0
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

        plt.title("Sphere Trajectory - Euler's Method")
        plt.xlabel("Horizontal Distance (m)")
        plt.ylabel("Height (m)")
        plt.grid()
        plt.show()

        # Graph 2: Velocity
        plt.figure(figsize=(10, 6))

        plt.plot(t, v[:, 0], label="Horizontal Velocity")
        plt.plot(t, v[:, 1], label="Vertical Velocity")

        plt.title("Sphere Velocity vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity (m/s)")
        plt.legend()
        plt.grid()
        plt.show()

    # RETURN RESULTS TABLE

    return table


# RUN SPHERE SIMULATION

sphere_trajectory = sphere(
    p0=(0, 0),
    v0=(10, 10),
    mass=3,
    radius=0.5,
    h=0.01
)

print(sphere_trajectory)

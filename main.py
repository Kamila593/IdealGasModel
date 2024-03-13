import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML


def arctan(y, x):
    if x < 0:
        return np.pi + np.arctan(y / x)
    elif x == 0:
        return 0
    else:
        if y < 0:
            return 2 * np.pi + np.arctan(y / x)
        else:
            return np.arctan(y / x)


class Box:

    def __init__(self, box_size=10, num_particles=100, time_step=0.1):
        super().__init__()
        self.size = box_size

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.plot_box_edges()

        self.ax.set_xlim(0, 3*box_size)
        self.ax.set_ylim(0, 3*box_size)
        self.ax.set_zlim(0, 3*box_size)

        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        self.particles = Particles(box_size, num_particles, time_step)

        # Create a slider for box size
        ax_box_size = plt.axes([0.1, 0.01, 0.65, 0.03])
        self.box_size_slider = Slider(ax_box_size, 'Box Size', 1, 30, valinit=box_size)
        self.box_size_slider.on_changed(self.update_box_edges)

    def plot_box_edges(self):
        # Plot the box edges
        box_edges = [
                      [0, 0, 0],
                      [0, 0, self.size],
                      [0, self.size, self.size],
                      [0, self.size, 0],
                      [0, 0, 0],
                      [self.size, 0, 0],
                      [self.size, 0, self.size],
                      [self.size, self.size, self.size],
                      [self.size, self.size, 0],
                      [self.size, 0, 0],
                      [0, 0, 0],
                      [0, self.size, 0],
                      [self.size, self.size, 0],
                      [self.size, 0, 0],
                      [self.size, 0, self.size],
                      [0, 0, self.size],
                      [0, self.size, self.size],
                      [self.size, self.size, self.size],
                      [self.size, 0, self.size]
                     ]

        box_edges = [[edge[i] for edge in box_edges] for i in range(3)]
        self.box_lines = self.ax.plot(*box_edges, color='black')

    def update_box_edges(self, new_size):
        self.size = new_size
        self.particles.update_box_size(new_size)
        for line in self.box_lines:
            line.remove()
        self.plot_box_edges()
        self.fig.canvas.draw_idle()

    def animate(self, frame_num):
        self.particles.move()
        self.scatter._offsets3d = (self.particles.positions[:, 0],
                                   self.particles.positions[:, 1],
                                   self.particles.positions[:, 2])
        return self.scatter,


class Particles:
    def __init__(self, box_size=10, num_particles=100, time_step=0.1):
        self.box_size = box_size
        self.num_particles = num_particles
        self.time_step = time_step
        self.positions = np.random.rand(num_particles, 3) * box_size    # Generate particles in random positions
        self.velocities = self.maxwell_boltzmann_velocity(num_particles, T, m)

    def maxwell_boltzmann_distribution(self, v, m, k, T):
        f = 4 * np.pi * (m / (2 * np.pi * k * T)) ** (3 / 2) * v ** 2 * np.exp((-m * v ** 2) / (2 * k * T))
        return f

    def maxwell_boltzmann_velocity(self, num_particles, T, m):
        # Generate velocities following the Maxwell-Boltzmann distribution
        v_max = np.sqrt(2 * k * T / m)
        v = np.linspace(0, 3 * v_max, 1000)
        prob = self.maxwell_boltzmann_distribution(v, m, k, T)
        # Sample velocities from the distribution
        velocities = np.random.choice(v, size=(num_particles, 3), p=prob/np.sum(prob))
        return velocities

    def update_box_size(self, new_size):
        self.box_size = new_size
        self.positions[self.positions > new_size] = new_size  # Limit particles to new box size

    def move(self):
        self.positions += self.velocities * self.time_step
        self.check_collisions()
        self.apply_boundary_conditions()

    def apply_boundary_conditions(self):
        # Implement bouncing behavior when particles touch the edges
        self.positions[self.positions < 0] = 0
        self.positions[self.positions > self.box_size] = self.box_size
        self.velocities[(self.positions == 0) | (self.positions == self.box_size)] *= -1

    def check_collisions(self):
        # Perform elastic collisions between particles
        for i in range(self.num_particles - 1):
            for j in range(i + 1, self.num_particles):

                distance = np.linalg.norm(self.positions[i] - self.positions[j])
                if distance <= 0.2:  # Check if particles collide

                    # Velocity vectors
                    v1 = self.velocities[i]
                    v2 = self.velocities[j]
                    r1 = np.linalg.norm(v1)
                    r2 = np.linalg.norm(v2)

                    # Move the beginning of coordinate system to the center of first particle
                    p1 = self.positions[i] - self.positions[i]
                    p2 = self.positions[j] - self.positions[i]

                    # Calculate vector between centers, which will be identified as rotation vector
                    k_vector = p2 - p1
                    k = np.linalg.norm(k_vector)
                    alpha = arctan(k_vector[1], k_vector[0])
                    beta = np.arccos(k_vector[2] / k)

                    # Calculate angles of velocity vectors
                    phi1 = arctan(v1[1], v1[0])
                    phi2 = arctan(v2[1], v2[0])
                    theta1 = np.arccos(v1[2] / r1)
                    theta2 = np.arccos(v2[2] / r2)

                    # Calculate velocities in rotated coordinate system
                    v1_r = np.array([
                        r1 * np.sin(theta1 + (np.pi / 2 - beta)) * np.cos(phi1 - alpha),
                        r1 * np.sin(theta1 + (np.pi / 2 - beta)) * np.sin(phi1 - alpha),
                        r1 * np.cos(theta1 + (np.pi / 2 - beta))
                    ])
                    v2_r = np.array([
                        r2 * np.sin(theta2 + (np.pi / 2 - beta)) * np.cos(phi2 - alpha),
                        r2 * np.sin(theta2 + (np.pi / 2 - beta)) * np.sin(phi2 - alpha),
                        r2 * np.cos(theta2 + (np.pi / 2 - beta))
                    ])

                    # Calculate normal components of velocities
                    v1_rn = v1_r[0]
                    v2_rn = v2_r[0]

                    # Calculate velocities after collision
                    v1__r = np.array([v2_rn, v1_r[1], v1_r[2]])
                    v2__r = np.array([v1_rn, v2_r[1], v2_r[2]])

                    # Calculate angles of velocity vectors after collision
                    r1__ = np.linalg.norm(v1__r)
                    r2__ = np.linalg.norm(v2__r)
                    phi1__ = arctan(v1__r[1], v1__r[0])
                    phi2__ = arctan(v2__r[1], v2__r[0])
                    theta1__ = np.arccos(v1__r[2] / r1__)
                    theta2__ = np.arccos(v2__r[2] / r2__)

                    # Transform velocities back to initial coordinate system
                    v1___ = np.array([
                        r1__ * np.sin(theta1__ - (np.pi / 2 - beta)) * np.cos(phi1__ + alpha),
                        r1__ * np.sin(theta1__ - (np.pi / 2 - beta)) * np.sin(phi1__ + alpha),
                        r1__ * np.cos(theta1__ - (np.pi / 2 - beta))
                    ])
                    v2___ = np.array([
                        r2__ * np.sin(theta2__ - (np.pi / 2 - beta)) * np.cos(phi2__ + alpha),
                        r2__ * np.sin(theta2__ - (np.pi / 2 - beta)) * np.sin(phi2__ + alpha),
                        r2__ * np.cos(theta2__ - (np.pi / 2 - beta))
                    ])

                    # Update velocities after collision
                    self.velocities[i] = v1___
                    self.velocities[j] = v2___


T = 80      # K
k = 1.38 * 10**(-23)        # J/K
m = 50 * 10**(-24)      # g

box = Box()
particles = Particles()

# Create the scatter plot of particles
box.scatter = box.ax.scatter(particles.positions[:, 0],
                             particles.positions[:, 1],
                             particles.positions[:, 2], color='red', s=20)

anim = FuncAnimation(box.fig, box.animate, frames=200, interval=5)
HTML(anim.to_jshtml())

plt.show()

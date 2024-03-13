## **Content**
1. [Program description](#program-description)
2. [Ideal Gas Model](#ideal-gas-model)
    * [Physical Chemistry Introduction](#physical-chemistry-introduction)
    * [Particles speeds](#particles-speeds)
    * [Elastic collisions](#elastic-collisions)
3. [How the program works](#how-the-program-works)
    * [Class *Box*](#class-box)
    * [Class *Particles*](#class-particles)
    * [Final part](#final-part)


## **Program description**

The following project is a computational model that visualizes the dynamics of an ideal gas. The model simulates particles enclosed inside a box, following Maxwell-Boltzmann distribution and elastic collision dynamics. User can modify the size of the box by moving slider displayed during runtime of the simulation, and observe the system from all angles thanks to dynamic view feature.



## **Ideal Gas**

### **Physical chemistry definition**

Ideal gas is a simplified model of the thermodynamic behaviour of gases. It is a set of identical particles that satisfies the following assumptions:
* particles are much smaller than their average distance apart, so they occupy negligible volume and can be considered to be point-like
* particles are in constant, chaotic motion
* particles don't interact with one another with intermolecular forces
* particles undergo elastic collisions between themselves and with the walls of the box
* average kinetic energy of the system is constant

Ideal gas model is useful for describing real gases, since it encapsulates principal gaseous properties. Therefore, in specific conditions real gases can be treated as ideal gases within reasonable tolerances.




### **Particle speeds**

According to the ideal gas model, individual molecules move at different speeds. Maxwell-Boltzmann distribution $f(v)$ is a probability density function that gives the probability of finding a particle with speed near $v$ and is given by the following equation:

$$
f(v) = \left(\frac{m}{2 \pi k T}\right)^\frac{3}{2} 4 \pi v^2 \exp{\left(-\frac{mv^2}{2kT}\right)}
$$

where:
* $m$ - particle mass
* $k$ - Boltzmann constant, $k$ = $\mathrm{1.38 \cdot 10^{-23}}$ $\frac{J}{K}$
* $T$ - temperature

As the equation states, particle speeds depend on temperature and molecular mass. Although the Maxwell-Boltzmann distribution gives detailed insight into the statistics of the occurrence of specific speeds in a system, the relationship of these parameters to speed can be explained directly. Indeed, the average particle speed in a system increases along with the temperature. The average kinetic energy of molecules, which is proportional to its temperature, directly influences average particle speed. Likewise, within the same temperature, light molecules on average move faster than heavy molecules. This can be explained based on the equation defining kinetic energy $E_k$,  $E_k~\mathrm{= \frac{1}{2}mv^2}$. If the molecules have the same average kinetic energy (same temperature), then lighter molecules must have a higher speed to compensate for their lower mass.




### **Elastic collisions**

Speeds are not fixed to individual particles, as they are constantly modified by collisions with other particles. A moving particle can be accelerated to very high speeds or slowed down by a collision with another particle. There are also collisions with walls of the box. They don't affect the magnitude of velocity vectors if the collision is elastic, but they change the direction of the vector.

In the ideal gas model all collisions are considered to be elastic. They occur when there's no loss in kinetic energy in the system as a result of the collision (no energy is scattered as heat). If a particle collides with a wall, the particle's velocity vector changes direction - the velocity component perpendicular to the wall is reversed, the rest of the components remain unchanged. If two particles collide with themselves, the situation is more complex. One way to calculate post-impact velocity vectors, which was used in this project, is explained below.

The center of the coordinate system is placed at the center of one of the molecules, then the coordinate system is rotated so that the x-axis passes through the center of the other molecule. This way, the collision can be reduced to one dimension and vector calculation can be performed on the vectors component parallel to the x-axis. Based on the components, the final velocity vectors after the collision are calculated and returned to the initial coordinate system.




## **How the program works**

After running the program, the simulation of ideal gas is displayed immediately. The user can rotate the box in every direction and modify the size of the box along with its volume. Gas enclosed inside the box expands with the increasing volume of the box, and contracts with decreasing volume of the box.

Code was written using OOP. Key features stating the behaviour of the box and particles were explained within two classes: *Box* and *Particles*.


### **Class** ***Box***

This class creates a container within which the simulated particles move. In the initialization, the introductory state of a box is defined. The method takes parameters such us box size, number of particles, and time step. The method sets up the initial size of the box, sets up a 3D plot for visualization, sets the limits for the plot axes based on the box size, and creates labels for the axes. It also creates a slider to dynamically adjust the size of the box during runtime. Additionally, it initializes a Particles object to represent the gas particles within the box. The following methods are defined within the class *Box*:
* *plot_box_edges* - method generates the coordinates for the edges of the box and plots them on the 3D plot; the box edges are defined as a series of vertices in 3D space with given length, forming the outline of the box
* *update_box_edges* - updates the size of the box when it's adjusted by a slider, redraws the box edges accordingly and changes the volume occupied by gas particles
* *animate* - is responsible for animating the movement of the gas particles within the box by updating their positions




### **Class** ***Particles***

At the beginning of a simulation, random positions of the particles are drawn. Although the particles are considered to be points, they are visualized as balls, so they are more visible. The number of particles is set to 100, but can be changed. All the particles get assigned velocities that are drawn based on the probability stated in the Maxwell-Boltzmann distribution. Program calculates new positions per each time step based on the previous position of a particle and based on the velocity. Particles move within a box and they can collide with sides of the box and with each other.

The collision definition must take into account the fact that time in the program is not continuous - the program calculates the positions of particles each selected time interval. Hence, the moment at which the particles come into contact with each other (a distance equal to two particle radii) cannot be considered a collision condition. This way, in most cases, in one time frame the particles would approach each other, and in the next time frame they would pass each other without colliding and move away. The collision condition must therefore define a distance greater than two particle radii. However, this method is not free from potential errors. If the collision condition is set to a distance that is relatively too large for a given time frame, then in the post-collision frame the program will 'think' that the particles are still within the collision range. This will trap the particles in a loop where they will bounce back and forth alternately.

Particles' movement is governed by the following methods:
* *maxwell_boltzmann_distribution* - calculates the probability density function for the Maxwell-Boltzmann distribution
* *maxwell_boltzmann_velocity* - generates velocities for the particles following the Maxwell-Boltzmann distribution
* *update_box_size* - adjusts volume of the gas to box limits - changing size of the box affects volume of the gas
* *apply_boundary_conditions* - is a function that adjusts velocities after a collision with a wall
* *check_collisions* - detects collisions between pairs of particles by calculating the distance between them; if particles are within a certain distance, they are considered to have collided, and their velocities are adjusted according to the laws of elastic collisions
* *move* - function updates the position of particles within time step of 0.1 s based on velocity of a given particle, and functions *apply_boundary_conditions* - and *check_collisions*



### **Final part**

At the end of the code, the values of temperature, Boltzmann's constant and particle's mass are defined. Values of temperature and particle's mass were adjusted in such a way that makes the movement during the simulation convenient to observe. The chosen initial values are: $T$ = 80 and $m$ = 1. However, it is possible to change the values of temperature and molecular mass. The way changing these parameters will affect particle speeds was explained in the *Particle speeds* section.

Objects 'box' and 'particles' are created based on classes 'Box' and 'Particles'. Having done that, everything is prepared for the simulation. The animation is displayed using HTML object.



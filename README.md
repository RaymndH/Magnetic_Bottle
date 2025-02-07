# Magnetic_Bottle
This was a project intended to incorporate a couple of techniques in numerical integration and 
computational physics in order to simulate the trapping effect of the fields of a magnetic mirror
on an arbitrary charged particle.

At this time, it consists of two files: Old Magnetic Bottle.ipynb and magnetic_bottle.py.
The former was based on the original project writeup for a class taken in the spring, 
the latter is revised code which amended several critical errors and has been streamlined. 

![new_animation](https://github.com/user-attachments/assets/c9efb783-4f47-46bf-bf07-3fc5147dbacc)

## Leapfrog algorithm
There are many different schemes for the numerical integration of forces on particles. One that is commonly used when
working with charged particles is the so-called Leapfrog algorithm. 

The naive method for integrating the equations of motion involves 2 variables for each particle, one storing the 3 dimensional 
position values and one storing the 3 dimensional velocity values. Since a charged particle in a field is subject to forces
and accelerations, one could, at each timestep, update the position based on the velocity times dt, and then the velocities as 
the accleration times dt. However, this is subject to instabilities. The leapfrog algorithm circumvents this by splitting the 
timesteps: dT is halved, and positions and velocities update at every other timestep.

This increases stability.

Structure of the main code: MPM.m

March 2018
Implement the TLS for MPM

INPUT DATA
1. MAterial Properties
2. Structure grid generation
3. Time variables
4. Compute Boundary nodes
5. Generate Particle varaibles (position and variables)
6. Generate Nodal variables
7. Initial condition

SOLVER (loop t<ft)
1. Reset nodal value to zero
2. Compute_Interpolator
3. Interpolate_Particle_To_Grid
4. Update momentum and boundary
5. Update_Particle Position and Velocity
6. Interpolate_velocity_back
7. Update_Stress (for paticles)


Additional for TLS:
1. Compute Gauss Location
2. Calculate the Integral
3. Compute the TLS fucntion for each cell
4. Gauss integration


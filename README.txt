Structure of the main code: MPM.m

March 2018
Implement the TLS for CPDI

INPUT DATA
1. Material Properties
2. Structure grid input and Grid generation
3. Time variables
4. Compute Boundary nodes
5. Generate Particle variables
6. Initial condition

SOLVER (loop t<ft)
1. Compute_Interpolator
2. Interpolate_Particle_To_Grid
3. Update momentum and boundary
4. Update_Particle Position and Velocity
5. Interpolate_velocity_back
6. Update_Stress (for paticles)


Additional for TLS for Interpolate_Particle_To_Grid and Interpolate_velocity_back:
1. Compute Gauss Location
2. Calculate the Integral
3. Compute the TLS fucntion for each cell
4. Gauss integration


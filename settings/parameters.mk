# --------------------------
# Parameters of the program
# --------------------------
# Units:
# cm      (centimeters)
# s       (seconds)

#Mesh Parameters:
Lx = 20                 	  #Size_x(cm)
Ly = 20                  	  #Size_y(cm)

#Simulation parameters:
DT = 5         #Dt(s)
T_FI = 100           #Final_time(s)
VIS = 0            #Visualization_steps

#FE parameters:
REF      = 4           #Refinements
ORDER    = 1           #Order
ABST_EX  = 0           #Absolute_Tolerance(Explicit)
RELT_EX  = 1E-8        #Relative_Tolerance(Explicit)
ITER_EX  = 200         #Iterations(Explicit)
ABST_IM  = 0           #Absolute_Tolerance(Implicit)
RELT_IM  = 1E-8        #Relative_Tolerance(Implicit)
ITER_IM  = 200         #Iterations(Implicit)
ABST_VEL = 0           #Absolute_Tolerance(Velocity)
RELT_VEL = 1E-8        #Relative_Tolerance(Velocity)
ITER_VEL = 200         #Iterations(Velocity)
ABST_SUN = 1E-5        #Absolute_Tolerance(SUNDIALS)
RELT_SUN = 1E-5        #Relative_Tolerance(SUNDIALS)
EPSILON  = 1E-9        #nEpsilon

#Vortex conditions:
Rx = 10				   #Position_x(cm)
Ry = 10				   #Position_y(cm)
STD = 2                #Standard_deviation(cm)
INT = 10			   #Intensity(cm²/s)
VISC = 0.15			   #Viscosity(cm²/s)

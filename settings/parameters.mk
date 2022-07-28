# --------------------------
# Parameters of the program
# --------------------------
# Units:
# cm      (centimeters)
# s       (seconds)

#Mesh Parameters:
Lx = 10                 	  #Size_x(cm)
Ly = 40                  	  #Size_y(cm)

#Simulation parameters:
DT = 0.001         #Dt(s)
T_FI = 10           #Final_time(s)
VIS = 100            #Visualization_steps

#FE parameters:
REF      = 4           #Refinements
ORDER    = 1           #Order
ABST_S   = 0           #Absolute_Tolerance(Solver)
RELT_S   = 1E-8        #Relative_Tolerance(Solver)
ITER_S   = 200         #Iterations(Solver)
ABST_SUN = 1E-5        #Absolute_Tolerance(SUNDIALS)
RELT_SUN = 1E-5        #Relative_Tolerance(SUNDIALS)
EPSILON  = 1E-9        #nEpsilon

#Vortex conditions:
Rx = 5				   #Position_x(cm)
Ry = 5				   #Position_y(cm)
STD = 2                #Standard_deviation(cm)
INT = 10			   #Intensity(cm²/s)
VISC = 0.15			   #Viscosity(cm²/s)

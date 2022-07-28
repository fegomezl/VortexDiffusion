#For Spack instalations
MFEM_INSTALL_DIR = /home/tapemess/Repositories/spack/opt/spack/linux-arch-sandybridge/gcc-12.1.0/mfem-develop-dc5wght7ajfaqp5se2zuqnhh6i7qob4u

#For source instalations
#MFEM_INSTALL_DIR = /opt/MFEM/mfem/build
#PETSC_INC = -I$(MFEM_INSTALL_DIR)/../../petsc/include
#SUNDIALS_INC = -I$(MFEM_INSTALL_DIR)/../../sundials/install/include

#If gmsh installed from source, else leave empty
GMSH_INSTALL = 

SHARE_DIR = /media/sf_Arch_data/
PARAVIEW_PATH =
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)

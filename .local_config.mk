#For Spack instalations
#MFEM_INSTALL_DIR = 

#For source instalations
#MFEM_INSTALL_DIR = /opt/MFEM/mfem/build
#PETSC_INC = -I$(MFEM_INSTALL_DIR)/../../petsc/include
#SUNDIALS_INC = -I$(MFEM_INSTALL_DIR)/../../sundials/install/include

#Local varibles
PROCCESORS = 1
SHARE_DIR = NULL

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)

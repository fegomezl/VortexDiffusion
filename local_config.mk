#For Spack instalations
MFEM_INSTALL_DIR = /home/boltzmann/.spack/opt/spack/linux-arch-westmere/gcc-12.2.0/mfem-develop-ot7cjjc7fwf373bsembbdjayqejcrgsd

#For source instalations
#MFEM_INSTALL_DIR = /opt/MFEM/mfem/build
#PETSC_INC = -I$(MFEM_INSTALL_DIR)/../../petsc/include
#SUNDIALS_INC = -I$(MFEM_INSTALL_DIR)/../../sundials/install/include

#Local varibles
PROCCESORS = 4
SHARE_DIR = /media/sf_ArchData/

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)

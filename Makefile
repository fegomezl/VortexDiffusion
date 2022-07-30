include local_config.mk
include settings/parameters.mk

#Compiling parameters
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS) $(PETSC_INC) $(SUNDIALS_INC)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main graph clean oclean

all: main

main: main.x
	@echo -e 'Running program ... \n'
	@$(RUN)$< -Lx $(Lx) -Ly $(Ly) \
		      -dt $(DT) -t_f $(T_FI) -v_s $(VIS) \
			  -ref $(REF) -o $(ORDER) -abstol_s $(ABST_S) -reltol_s $(RELT_S) -iter_s $(ITER_S) \
			  -abstol_sun $(ABST_SUN) -reltol_sun $(RELT_SUN) -eps $(EPSILON) \
			  -Rx $(Rx) -Ry $(Ry) -sigma $(STD) -gamma $(INT) -visc $(VISC)
	@echo -e '\n\nDone!\n'

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/vortex_diffusion
	@cp -r results/graph $(SHARE_DIR)/vortex_diffusion
	@echo -e 'Done!'
endif

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x results/graph/*

oclean:
	@rm -rf .objects/*.o

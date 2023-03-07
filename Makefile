FC = gfortran
CC = gcc
LINALG = -llapack -lblas
MODCMD = -J

OBJDIR = obj
SRCDIR = src
MODDIR = mod

MAINDIR=main
SAODIR=system_and_options
TYPDIR=types
INTDIR=interfaces
TASKDIR=tasks
PHYSDIR=physics
CALCDIR=calculations
PAWDIR=libpaw


MPI_COMPILE_FLAGS = $(shell mpif90 --showme:compile)
MPI_LINK_FLAGS = $(shell mpif90 --showme:link)
MPI_INC_DIR= $(shell mpif90 --showme:incdirs)

#FFTW Install with mpi, threads enabled
FFTW_DIR=/usr/local/opt/fftw/
FFTW_COMPILE_FLAGS =-I$(FFTW_DIR)include -I$(FFTW_DIR)lib
FFTW_LINK_FLAGS = -L$(FFTW_DIR)lib -lfftw3_mpi -lfftw3 -lfftw3_threads -lm

XC_DIR=/usr/local/opt/libxc/
XC_COMPILE_FLAGS =-I$(XC_DIR)include -I$(XC_DIR)lib
XC_LINK_FLAGS = -L$(XC_DIR)lib -lxcf03 -lxc

LA_DIR=/usr/local/opt/lapack/
LA_COMPILE_FLAGS =-I$(LA_DIR)include -I$(LA_DIR)lib -llapack -lblas

FFT=fftw
parallel=mpi
xc=libxc
linear_algebra=lapack

INC= -I$(MODDIR)/$(PAWDIR) -I$(MODDIR)/$(TYPDIR) -I$(MODDIR)/$(INTDIR) -I$(MODDIR)/$(SAODIR) -I$(MODDIR)/$(TASKDIR) -I$(MODDIR)/$(PHYSDIR) -I$(MODDIR)/$(CALCDIR) -I$(MODDIR)/$(MAINDIR) $(MPI_COMPILE_FLAGS) $(FFTW_COMPILE_FLAGS) $(XC_COMPILE_FLAGS) 

MOD_DIRECTORIES= $(MODDIR)/$(MAINDIR) $(MODDIR)/$(CALCDIR) $(MODDIR)/$(PHYSDIR) $(MODDIR)/$(TASKDIR) $(MODDIR)/$(SAODIR) $(MODDIR)/$(INTDIR) $(MODDIR)/$(TYPDIR) $(MODDIR)/$(PAWDIR)
OBJ_DIRECTORIES= $(OBJDIR)/$(MAINDIR) $(OBJDIR)/$(CALCDIR) $(OBJDIR)/$(PHYSDIR) $(OBJDIR)/$(TASKDIR) $(OBJDIR)/$(SAODIR) $(OBJDIR)/$(INTDIR) $(OBJDIR)/$(TYPDIR) $(OBJDIR)/$(PAWDIR)

LINK = $(LINALG) $(MPI_LINK_FLAGS) $(XC_LINK_FLAGS) $(FFTW_LINK_FLAGS)

$(shell   mkdir -p $(MOD_DIRECTORIES)) 
$(shell   mkdir -p $(OBJ_DIRECTORIES)) 

OBJPAW= \
		$(OBJDIR)/$(PAWDIR)/libpaw_libxc.o\
		$(OBJDIR)/$(PAWDIR)/m_libpaw_defs.o\
		$(OBJDIR)/$(PAWDIR)/m_libpaw_mpi.o\
		$(OBJDIR)/$(PAWDIR)/m_libpaw_tools.o\
		$(OBJDIR)/$(PAWDIR)/m_libpaw_libxc.o\
		$(OBJDIR)/$(PAWDIR)/m_paral_atom.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_numeric.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_sphharm.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_io.o\
		$(OBJDIR)/$(PAWDIR)/m_pawang.o\
		$(OBJDIR)/$(PAWDIR)/m_pawrad.o\
		$(OBJDIR)/$(PAWDIR)/m_pawtab.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_an.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_ij.o\
		$(OBJDIR)/$(PAWDIR)/m_pawfgrtab.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_finegrid.o\
		$(OBJDIR)/$(PAWDIR)/m_pawcprj.o\
		$(OBJDIR)/$(PAWDIR)/m_pawrhoij.o\
		$(OBJDIR)/$(PAWDIR)/m_pawdij.o\
		$(OBJDIR)/$(PAWDIR)/m_pawxc.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_atom.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_gaussfit.o\
		$(OBJDIR)/$(PAWDIR)/m_pawxmlps.o\
		$(OBJDIR)/$(PAWDIR)/m_pawpsp.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_init.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_occupancies.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_correlations.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_denpot.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_nhat.o\
		$(OBJDIR)/$(PAWDIR)/m_paw_onsite.o


OBJTYPE= \
		$(OBJDIR)/$(TYPDIR)/types.o\
		$(OBJDIR)/$(TYPDIR)/parallel_type_$(parallel).o\
		$(OBJDIR)/$(TYPDIR)/shared_type.o\
		$(OBJDIR)/$(TYPDIR)/constants.o\
		$(OBJDIR)/$(TYPDIR)/numerics.o\
		$(OBJDIR)/$(TYPDIR)/odp_type_$(FFT).o\
		$(OBJDIR)/$(TYPDIR)/element_type.o\
		$(OBJDIR)/$(TYPDIR)/atom_type.o\
		$(OBJDIR)/$(TYPDIR)/grids_type_$(FFT).o\
		$(OBJDIR)/$(TYPDIR)/system_type.o\
		$(OBJDIR)/$(TYPDIR)/xc_type_$(xc).o\
		$(OBJDIR)/$(TYPDIR)/tf_type.o\
		$(OBJDIR)/$(TYPDIR)/td_type.o\
		$(OBJDIR)/$(TYPDIR)/ewald_type.o\
		$(OBJDIR)/$(TYPDIR)/KG_type.o\
		$(OBJDIR)/$(TYPDIR)/DOS_type.o\
		$(OBJDIR)/$(TYPDIR)/tasks_type.o\
		$(OBJDIR)/$(TYPDIR)/stochastic_type.o\
		$(OBJDIR)/$(TYPDIR)/stopping_type.o\
		$(OBJDIR)/$(TYPDIR)/simulation_type.o

OBJINT = \
		$(OBJDIR)/$(INTDIR)/utils.o\
		$(OBJDIR)/$(INTDIR)/$(FFT).o\
		$(OBJDIR)/$(INTDIR)/xc_interface_$(xc).o\
		$(OBJDIR)/$(INTDIR)/parallel_$(parallel).o\
		$(OBJDIR)/$(INTDIR)/lapack.o\
		$(OBJDIR)/$(INTDIR)/linalg.o\
		$(OBJDIR)/$(INTDIR)/splines.o

OBJSAO = \
        $(OBJDIR)/$(SAODIR)/inputs.o\
		$(OBJDIR)/$(SAODIR)/grids_$(FFT).o\
		$(OBJDIR)/$(SAODIR)/odp_$(FFT).o

OBJTASK = \
		$(OBJDIR)/$(TASKDIR)/operations_3D.o\
		$(OBJDIR)/$(TASKDIR)/Spherical_Harmonics.o\
		$(OBJDIR)/$(TASKDIR)/field_gradient.o\
		$(OBJDIR)/$(TASKDIR)/distribute_field.o\
		$(OBJDIR)/$(TASKDIR)/FFT2Matrix.o\
		$(OBJDIR)/$(TASKDIR)/Fourier.o\
		$(OBJDIR)/$(TASKDIR)/read_in.o

OBJPHYS = \
		$(OBJDIR)/$(PHYSDIR)/Ewald.o\
		$(OBJDIR)/$(PHYSDIR)/Hartree.o\
		$(OBJDIR)/$(PHYSDIR)/Kinetic.o\
		$(OBJDIR)/$(PHYSDIR)/External_Field.o\
		$(OBJDIR)/$(PHYSDIR)/XC.o\
		$(OBJDIR)/$(PHYSDIR)/Density.o\
		$(OBJDIR)/$(PHYSDIR)/HGH_Pseudopotentials.o\
		$(OBJDIR)/$(PHYSDIR)/PAW_Pseudopotentials.o\
		$(OBJDIR)/$(PHYSDIR)/Local_Ion.o\
		$(OBJDIR)/$(PHYSDIR)/Non_Local_Ion.o\
		$(OBJDIR)/$(PHYSDIR)/Apply_Hamiltonian.o\
		$(OBJDIR)/$(PHYSDIR)/Non_Local_Ion_EnFoSt.o\
		$(OBJDIR)/$(PHYSDIR)/Chebychev.o\
		$(OBJDIR)/$(PHYSDIR)/Stochastic.o\
		$(OBJDIR)/$(PHYSDIR)/Thermostat.o\
		$(OBJDIR)/$(PHYSDIR)/Thomas_Fermi.o\
		$(OBJDIR)/$(PHYSDIR)/Entropy.o\
		$(OBJDIR)/$(PHYSDIR)/Chemical_Potential.o

OBJCALC = \
		$(OBJDIR)/$(CALCDIR)/allocate_memory.o\
		$(OBJDIR)/$(CALCDIR)/Orbital_Free_Min.o\
		$(OBJDIR)/$(CALCDIR)/Update_Potentials_Energies.o\
		$(OBJDIR)/$(CALCDIR)/initialization_finalization.o\
		$(OBJDIR)/$(CALCDIR)/Eigen_LOBPCG.o\
		$(OBJDIR)/$(CALCDIR)/Eigen_ChebFilter.o\
		$(OBJDIR)/$(CALCDIR)/Eigen_RMM_DIIS.o\
		$(OBJDIR)/$(CALCDIR)/Stress_Pressure.o\
		$(OBJDIR)/$(CALCDIR)/SCF.o\
		$(OBJDIR)/$(CALCDIR)/Stopping_Power.o\
		$(OBJDIR)/$(CALCDIR)/Current.o\
		$(OBJDIR)/$(CALCDIR)/Time_Dependent.o\
		$(OBJDIR)/$(CALCDIR)/Kubo_Greenwood.o\
		$(OBJDIR)/$(CALCDIR)/DensityofStates.o\
		$(OBJDIR)/$(CALCDIR)/Simple_Tests.o

OBJMAIN = \
		$(OBJDIR)/$(MAINDIR)/main.o 

$(OBJDIR)/$(PAWDIR)/%.o: $(SRCDIR)/$(PAWDIR)/%.c
	$(CC) $(INC) $(CFLAG) -o $@ -c $< 

$(OBJDIR)/$(PAWDIR)/%.o: $(SRCDIR)/$(PAWDIR)/%.F90
	$(FC) $(INC) $(FFLAG) -o $@ -c $<  $(MODCMD)$(MODDIR)/$(PAWDIR)/

##Compilation order types > interfaces(libraries) > systems & option > tasks > physics > main
$(OBJDIR)/$(TYPDIR)/%.o: $(SRCDIR)/$(TYPDIR)/%.f90
	$(FC) $(INC)  $(FFLAG)  -o $@ -c $<  $(MODCMD)$(MODDIR)/$(TYPDIR)/

$(OBJDIR)/$(INTDIR)/%.o: $(SRCDIR)/$(INTDIR)/%.f90
	$(FC) $(INC)  $(FFLAG) -o $@ -c $<  $(MODCMD)$(MODDIR)/$(INTDIR)/

$(OBJDIR)/$(SAODIR)/%.o: $(SRCDIR)/$(SAODIR)/%.f90
	$(FC) $(INC)  $(FFLAG) -o $@ -c $<  $(MODCMD)$(MODDIR)/$(SAODIR)/

$(OBJDIR)/$(TASKDIR)/%.o: $(SRCDIR)/$(TASKDIR)/%.f90
	$(FC) $(INC)  $(FFLAG) -o $@ -c $<  $(MODCMD)$(MODDIR)/$(TASKDIR)/

$(OBJDIR)/$(PHYSDIR)/%.o: $(SRCDIR)/$(PHYSDIR)/%.f90
	$(FC) $(INC)  $(FFLAG) -o $@ -c $<  $(MODCMD)$(MODDIR)/$(PHYSDIR)/

$(OBJDIR)/$(CALCDIR)/%.o: $(SRCDIR)/$(CALCDIR)/%.f90
	$(FC) $(INC)  $(FFLAG) -o $@ -c $<  $(MODCMD)$(MODDIR)/$(CALCDIR)/

$(OBJDIR)/$(MAINDIR)/%.o: $(SRCDIR)/$(MAINDIR)/main.f90
	$(FC) $(INC)  $(FFLAG) -o $@ -c $< 

gnu:  FC = gfortran 
gnu:  CC = gcc
gnu:  FFLAG =  -O2 -mcmodel=medium -std=f2018
gnu:  CFLAG= -O2
gnu:  LDFLAGS = $(FFLAG)
gnu:  MODCMD = -J
gnu:  main.exe

gnu_mkl: FC = gfortran
gnu_mkl: CC = /usr/bin/gcc
gnu_mkl: LINALG =  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
gnu_mkl: FFLAG= -O3 -I${MKLROOT}/include -mcmodel=medium
gnu_mkl: CFLAG= -O3 -I${MKLROOT}/include -DMKL_LP64
gnu_mkl: MODCMD = -J
gnu_mkl: main.exe

gnu_debug:  FC = gfortran
gnu_debug:  CC = /usr/bin/gcc
gnu_debug:  FFLAG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow 
gnu_debug:  LDFLAGS = $(FFLAG)
gnu_debug:  CFLAG= -g
gnu_debug: MODCMD = -J
gnu_debug:  main.exe

main.exe:  $(OBJPAW) $(OBJTYPE) $(OBJINT) $(OBJSAO) $(OBJTASK) $(OBJPHYS) $(OBJCALC) $(OBJMAIN)
			$(FC) $(LDFLAGS) $(INC) -o main.exe $(OBJPAW) $(OBJTYPE) $(OBJINT) $(OBJSAO) $(OBJTASK) $(OBJPHYS) $(OBJCALC) $(OBJMAIN) $(LINK)

	#Copy module files into src directories for visual studios autodebugging
	$(shell   cp $(MODDIR)/$(PAWDIR)/*.mod $(SRCDIR)/$(PAWDIR)/ )
	$(shell   cp $(MODDIR)/$(TYPDIR)/*.mod $(SRCDIR)/$(TYPDIR)/ )
	$(shell   cp $(MODDIR)/$(INTDIR)/*.mod $(SRCDIR)/$(INTDIR)/ )
	$(shell   cp $(MODDIR)/$(SAODIR)/*.mod $(SRCDIR)/$(SAODIR)/ )
	$(shell   cp $(MODDIR)/$(TASKDIR)/*.mod $(SRCDIR)/$(TASKDIR)/ )
	$(shell   cp $(MODDIR)/$(PHYSDIR)/*.mod $(SRCDIR)/$(PHYSDIR)/ )
	$(shell   cp $(MODDIR)/$(CALCDIR)/*.mod $(SRCDIR)/$(CALCDIR)/ )
	$(shell   cp $(MPI_INC_DIR)/mpi.mod $(SRCDIR)/$(INTDIR)/ )

clean :
	rm -rf $(OBJDIR) $(MODDIR) $(LIBDIR) $(SRCDIR)/$(PAWDIR)/*.mod $(SRCDIR)/$(TYPDIR)/*.mod $(SRCDIR)/$(INTDIR)/*.mod $(SRCDIR)/$(SAODIR)/*.mod $(SRCDIR)/$(TASKDIR)/*.mod $(SRCDIR)/$(PHYSDIR)/*.mod $(SRCDIR)/$(CALCDIR)/*.mod 
	rm main.exe

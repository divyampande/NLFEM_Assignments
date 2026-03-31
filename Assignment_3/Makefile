# ==========================================
# Fortran FEA Solver Makefile
# ==========================================

# 1. Variables (The Factory Settings)
FC = gfortran
FFLAGS = -O3 -Wall
TARGET = fea_solver

# 2. Source Files (ORDER MATTERS HERE!)
SRCS = kinds_mod.f90 material_mod.f90 solver_mod.f90 main.f90

# 3. Object Files (Swap .f90 for .o)
OBJS = $(SRCS:.f90=.o)

# ==========================================
# Build Rules
# ==========================================

# Default rule: Build the final executable
all: $(TARGET)

# Link the object files to create the executable
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) -llapack -lblas

# Rule to compile individual .f90 files into .o files
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Rule to clean up the factory floor (delete .o, .mod, and the executable)
clean:
	del /Q /F *.o *.mod $(TARGET).exe

# Rule to build AND run the code immediately
run: $(TARGET)
	python run_sim.py
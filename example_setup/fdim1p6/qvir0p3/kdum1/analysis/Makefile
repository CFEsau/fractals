# Compiler options
TARGET	    :=  reduce	
F90         :=  gfortran
F77	    :=  gfortran
FFLAGS      :=  -O3
VPATH       :=  src
BINDIR	    :=  bin
OBJDIR      :=  obj
LINKER 	    := gfortran -o

# Define file extensions
.SUFFIXES: .f90 .f .o .mod

# change these to set the proper directories where each files should be
OBJECTS  =  $(OBJDIR)/modules.o $(OBJDIR)/read_sl_out.o $(OBJDIR)/read_sl_sub.o $(OBJDIR)/reduce.o $(OBJDIR)/heapsort.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINKER) $@ $(OBJECTS)

$(OBJDIR)/%.o: $(VPATH)/%.f90 | $(OBJDIR)
	$(F90) $(FFLAGS) -c $^ -o $@ -J$(OBJDIR)

$(OBJDIR): FORCE
	mkdir -p $(OBJDIR)

FORCE:

clean:
	rm -rf $(OBJDIR)
	rm -rf $(TARGET)


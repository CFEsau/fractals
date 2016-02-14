# Compiler options
F90         :=  gfortran
F77	    :=  gfortran
FFLAGS      :=  -O3
VPATH       :=  src
OBJDIR      :=  obj

# Define file extensions
.SUFFIXES: .f90 .f .o .mod 

# All modules
OBJS        :=  $(OBJDIR)/modules.o $(OBJDIR)/plummer.o $(OBJDIR)/ran2.o $(OBJDIR)/assign_binary.o $(OBJDIR)/orientate.o $(OBJDIR)/makefractal.o  $(OBJDIR)/maschberger_imf.o

# Build rules
all: initials | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
	mkdir -p output

$(OBJDIR)/%.o: $(VPATH)/%.f90 | $(OBJDIR)
	$(F90) $(FFLAGS) -c $^ -o $@ -J$(OBJDIR)

initials.o: initials.f90 | $(OBJDIR)
	$(F90) $(FFLAGS) -I$(OBJDIR) -c $^ -o $@ -J$(OBJDIR)

initials: $(OBJS) $(OBJDIR)/initials.o
	$(F90) $(FFLAGS) -o initials $(OBJS) $(OBJDIR)/initials.o

clean:
	@rm -rf $(OBJDIR)
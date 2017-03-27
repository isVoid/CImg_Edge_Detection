# 'make depend' uses makedepend to automatically generate dependencies 
#
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

SWTYPE = CImg
SWNAME = HOUGH_LINE

# define the C compiler to use
CC = g++

# define any compile-time flags
SCFLAGS = -s
CFLAGS = -Wall -g

# define any directories containing header files other than /usr/include
#
INCLUDES = -I./include

LIBS = 

LDFLAGS = -lpthread -lX11 -ljpeg

# define the C source files
MAINSRCS = ./src/hough.cpp ./src/canny.h ./src/canny.cpp

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
MAINOBJS = $(MAINSRCS:.c=.o)

# define the executable file 
MAIN = $(SWNAME)

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN) $(AUTO)
		@echo  $(SWTYPE) software $(MAIN) has been compiled

$(MAIN): $(MAINOBJS) 
		$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(MAINOBJS) $(LIBS) $(LDFLAGS)

silent: $(OBJS) 
		$(CC) $(SCFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LIBS) $(LDFLAGS)
		@g++ note


# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
		$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
		$(RM) *.o *~ $(MAIN) $(AUTO)

depend: $(SRCS)
		makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

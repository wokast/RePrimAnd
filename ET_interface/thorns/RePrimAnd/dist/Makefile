# Simple Makefile for use with the ET's auto-buid mechanism
SUBDIRS:=$(wildcard */)
SOURCES:=$(wildcard $(patsubst %,%/*.cc,$(SUBDIRS)))
INCLUDES:=$(wildcard $(patsubst %,%/include/*.h,$(SUBDIRS)))
OBJECTS:=$(patsubst %.cc,%.o,$(SOURCES))

CXXFLAGS+=$(patsubst %,-I%/include,$(SUBDIRS))
CXXFLAGS+=$(patsubst %,-I%,$(BOOST_INC_DIRS) $(GSL_INC_DIRS) $(HDF5_INC_DIRS))

all: libRePrimAnd.a $(INCLUDES)

install: all
	mkdir -p $(DESTDIR)/include/reprimand $(DESTDIR)/lib
	cp libRePrimAnd.a $(DESTDIR)/lib
	cp $(INCLUDES) $(DESTDIR)/include/reprimand

libRePrimAnd.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

clean:
	$(RM) $(OBJECTS) libRePrimAnd.a

.PHONY: all install clean

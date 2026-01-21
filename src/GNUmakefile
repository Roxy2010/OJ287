.SUFFIXES: .cxx .d .o

# C++ compiler.
CXX = g++

# Compiler and linker options.
CXXFLAGS = -O3 -Wall
LINKFLAGS = 

# Compiler flags for dependency generation.
# Comment it out if your compiler cannot generate dependencies.
DEPFLAGS = -MM

LIBS = -lm
INCLUDE = -I. -I../include/

BINPROGS = cbwaves 
O = cbwaves.o 
D = cbwaves.d 

prefix=/usr/local/cbwaves

all: $(BINPROGS)

cbwaves: cbwaves.o
	$(CXX) $(LINKFLAGS) -o cbwaves cbwaves.o ConfigFile.o

clean::
	rm -f $(BINPROGS) $(O) $(D) ConfigFile.o

distclean:: clean
	rm -f $(O:.o=.d)

install:: cbwaves
	mkdir -p ${prefix}
	mkdir -p ${prefix}/bin
	mkdir -p ${prefix}/man
	mkdir -p ${prefix}/etc
	mkdir -p ${prefix}/doc
	mkdir -p ${prefix}/examples
	cp cbwaves ${prefix}/bin/
	chmod +x ${prefix}/bin/cbwaves
	cp ../man/cbwaves.1 ${prefix}/man
	cp ../etc/cbwaves.ini.tpl ${prefix}/etc
	cp ../examples/cbwgen.pl ${prefix}/examples/
	cp ../doc/README ${prefix}/doc

tgz:: distclean
	(cd .. ; rm -rf /tmp/cbwaves-`cat VERSION.txt`-`cat RELEASE`)
	(cd .. ; mkdir -p /tmp/cbwaves-`cat VERSION.txt`-`cat RELEASE`)
	(cd .. ; cp -R * /tmp/cbwaves-`cat VERSION.txt`-`cat RELEASE`/)
	(cd .. ; tar -C /tmp -czf /tmp/cbwaves-`cat VERSION.txt`-`cat RELEASE`.tgz cbwaves-`cat VERSION.txt`-`cat RELEASE`)
	(cd .. ; rm -rf /tmp/cbwaves-`cat VERSION.txt`-`cat RELEASE`)

rpm:: tgz
	(rm -rf /tmp/cbwavesrpm/)
	(mkdir -p /tmp/cbwavesrpm/SOURCES) 
	(mkdir -p /tmp/cbwavesrpm/BUILD) 
	(mkdir -p /tmp/cbwavesrpm/RPMS/x86_64) 
	(mkdir -p /tmp/cbwavesrpm/SRPMS) 
	(mkdir -p /tmp/cbwavesrpm/RPM/i386) 
	(cd .. ; cp /tmp/cbwaves-`cat VERSION.txt`-`cat RELEASE`.tgz /tmp/cbwavesrpm/SOURCES);
	(cd .. ; rpmbuild -v -ba cbwaves.spec)
	
.cxx.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c ConfigFile.cxx $<

.cxx.d:
	@echo "Creating dependency file $@"
	$(CXX) $(DEPFLAGS) $(INCLUDE) $< | sed 's|$*\.o: |$*\.o $*\.d: |' > $@        

-include $(O:.o=.d)

CPPFLAGS = -I.
CXXFLAGS = -Wall -g

OBJECTS = Tissue.o main.o
LIBS = fwk/BaseCollection.o fwk/BaseNotifiee.o fwk/Exception.o
LOG = infectionLog.txt

asgn1:	$(OBJECTS) $(LIBS)
	$(CXX) $(CXXFLAGS) -o asgn1 $(OBJECTS) $(LIBS)

clean:
	rm -f asgn1 $(OBJECTS) $(LIBS) *~ $(LOG)

Tissue.o: Tissue.cpp Tissue.h
main.o: main.cpp

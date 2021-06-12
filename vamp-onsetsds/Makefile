
CXXFLAGS	:= -I../vamp-plugin-sdk -O3 -Wall
#CXXFLAGS	:= -I../vamp-plugin-sdk -g -Wall -march=pentium4 -msse -msse2 -ffast-math
#CXXFLAGS	:= -I../vamp-plugin-sdk -O3 -Wall -march=pentium4 -msse -msse2 -fomit-frame-pointer -ffast-math

vamp-onsetsds.so:	onsetsdsplugin.o onsetsds/onsetsds.o
	g++ -shared $^ -o $@ -L../vamp-plugin-sdk/vamp-sdk -Wl,-Bstatic -lvamp-sdk -Wl,-Bdynamic -lpthread

clean:	
	rm -f *.o */*.o

onsetsdsplugin.o: onsetsdsplugin.h onsetsds/onsetsds.h 
onsetsds/onsetsds.o: onsetsds/onsetsds.h 

LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
 
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

canonicRecursion: canonicRecursion.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

clean:
	rm *.o canonicRecursion

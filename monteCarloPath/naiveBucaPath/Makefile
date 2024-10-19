LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
 
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

naiveBucaPath: naiveBucaPath.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

clean:
	rm *.o naiveBucaPath

LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
 
%.o : %.cpp
	g++ -Wall -c $< ${INCS}

idealHarmonic: idealHarmonic.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

clean:
	rm *.o idealHarmonic

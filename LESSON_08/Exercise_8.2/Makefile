LIBS:=`root-config --libs`
INCS:=`root-config --cflags`
 
%.o : %.cpp
	g++ -Wall -c $< 

qmc1d: qmc1d.o random.o 
	g++ -O3 -Wall -o $@ $^ 

clean:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat

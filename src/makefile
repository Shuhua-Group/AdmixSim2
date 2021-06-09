AdmixSim2: Segment.o Utils.o Haplotype.o Individual.o Generation.o Population.o Param.o AdmixSim2.o
	g++ Segment.o Utils.o Haplotype.o Individual.o Generation.o Population.o Param.o AdmixSim2.o -o AdmixSim2 -O3 

Segment.o: Segment.cpp Segment.hpp
	g++ -c Segment.cpp -O3

Utils.o: Utils.cpp Utils.hpp
	g++ -c Utils.cpp -O3

Haplotype.o: Haplotype.cpp Haplotype.hpp
	g++ -c Haplotype.cpp -O3

Individual.o: Individual.cpp Individual.hpp
	g++ -c Individual.cpp -O3

Generation.o: Generation.cpp Generation.hpp
	g++ -c Generation.cpp -O3

Population.o: Population.cpp Population.hpp
	g++ -c Population.cpp -O3

Param.o: Param.cpp Param.hpp
	g++ -c Param.cpp -O3

AdmixSim2.o: AdmixSim2.cpp
	g++ -c AdmixSim2.cpp -O3

clean:
	rm *.o AdmixSim2


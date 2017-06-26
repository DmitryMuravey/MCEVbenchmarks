TARGET = MCEVbenchmark

all: $(TARGET)

clean:
	rm -rf $(TARGET) *.o

$(TARGET).o:
	g++ -c -o $(TARGET).o $(TARGET).cc -g -O3 -mtune=native -std=c++11

$(TARGET): MCEVbenchmark.o
	g++ -o $(TARGET) $(TARGET).o -lgsl -lgslcblas


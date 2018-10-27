INCLUDES := -Iincludes
LIBS := -lGL -lglut
FLAGS := -fopenmp -g `pkg-config --cflags --libs opencv`
CC := g++
OUTPUT := pt_cpu
all: build

build:
	$(CC) main.cpp $(FLAGS) $(INCLUDES) $(LIBS) -o $(OUTPUT)
clean:
	rm -f $(OUTPUT)
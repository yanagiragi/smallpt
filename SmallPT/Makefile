INCLUDES := -Iincludes
LIBS := -lGL -lglut
FLAGS := -fopenmp -g `pkg-config --cflags --libs opencv`
CC := g++
all: build

build:
	$(CC) main.cpp $(FLAGS) $(INCLUDES) $(LIBS) -o pt_cpu

INCLUDES := -Iincludes
LIBS := -lGL -lglut
FLAGS := -O3 -fopenmp -g `pkg-config --cflags --libs opencv`

all: build

build:
	g++ main.cpp $(FLAGS) $(INCLUDES) $(LIBS) -o pt_cpu
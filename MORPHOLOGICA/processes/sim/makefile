PROGRAM = process
SRCS    = $(shell find . -name "*.cpp")
OBJS    = $(SRCS:.cpp=.o)
CC      = g++
CFLAGS  = -I/opt/X11/include -O2 -g `pkg-config opencv --cflags` -I/usr/include -DARMA_DONT_USE_WRAPPER
LDFLAGS = -L/usr/X11R6/lib -lX11 -lGL -framework OpenGL -llapack -lblas `pkg-config opencv --libs`


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o:%.cpp
	$(CC) -c -o $@ $(CFLAGS) $<

#.PHONY: clean dist-clean

clean:
	rm $(OBJS)

dist-clean: clean
	rm $(PROGRAM)

doc:
	doxygen doxygen.conf

doc-clean:
	rm -rf ./html









#g++ -o process process.cpp -I /opt/X11/include lib/display.cpp lib/brick.cpp ../lib/sockserve.cpp lib/body.cpp -L /opt/X11/lib -lX11 -lGL -framework OpenGL -framework GLUT `pkg-config opencv --libs --cflags` -I/usr/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

#LDFLAGS = -L/opt/local/lib -L/usr/X11R6/lib -L/opt/X11/lib -lX11 -lGL -framework OpenGL -llapack -lblas -lArmadillo -framework Accelerate `pkg-config opencv --libs`

CFLAGS=-w

all: _lm.so

_lm.so: 	lm.o lm_wrap.o
	g++ $(CFLAGS) -shared -o _lm.so lm.o lm_wrap.o -O1 -larmadillo -lcminpack

lm.o: 		lm.cpp
	g++ $(CFLAGS) -c -x c++ lm.cpp

lm_wrap.o: 	lm_wrap.c
	g++ -c lm_wrap.c -I/usr/include/python2.7/ -I/usr/local/lib/python2.7/dist-packages/numpy/core/include/

clean:
	rm lm_wrap.c *.o _lm.so *.pyc

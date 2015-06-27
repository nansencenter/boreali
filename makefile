CFLAGS=-I/usr/include/cminpack-1/ -L/usr/lib64/ -I/usr/include/python2.6/ -I/usr/lib64/python2.6/site-packages/numpy/core/include/ -w -fPIC

all: _lm.so

_lm.so: 	lm.o lm_wrap.o
	g++ $(CFLAGS) -shared -o _lm.so lm.o lm_wrap.o -L/opt/cminpack/1.3.0/ -O1 -lcminpack

lm.o: 		lm.cpp
	g++ $(CFLAGS) -c -x c++ lm.cpp -I/opt/cminpack/1.3.0/ -I/usr/include/cminpack-1/

lm_wrap.o: 	lm_wrap.c
	g++ $(CFLAGS) -c lm_wrap.c -I/usr/include/python2.7/ -I/usr/local/lib/python2.7/dist-packages/numpy/core/include/

clean:
	rm *.o _lm.so *.pyc

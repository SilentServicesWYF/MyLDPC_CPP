cc = g++
CFLAGS = -g

all: decoding

decoding: auxil.o decoding.o llrv.o gfcalu.o
	$(cc) $^ -o $@

decoding.o: decoding.cpp auxil.h llrv.h gfcalu.h
	$(cc) $(CFLAGS) -c $< -o $@

auxil.o: auxil.cpp auxil.h
	$(cc) $(CFLAGS) -c $< -o $@

llrv.o: llrv.cpp llrv.h gfcalu.h
	$(cc) $(CFLAGS) -c $< -o $@

gfcalu.o: gfcalu.cpp gfcalu.h
	$(cc) $(CFLAGS) -c $< -o $@

.PHONY:clean
clean:
	rm decoding decoding.o auxil.o llrv.o gfcalu.o
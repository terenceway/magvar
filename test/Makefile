.PHONY: all check clean

INCLUDE=..

CFLAGS=-I$(INCLUDE)

all : test_magvar

OBJS = test_magvar.o av_mag.o wmm.o

test_magvar : $(OBJS)
	$(CC) -o $@ $(OBJS)

wmm.o : ../wmm.c ../av_geo.h
	$(CC) $(CFLAGS) -c ../wmm.c

av_mag.o : ../av_mag.c ../av_mag.h ../av_geo.h
	$(CC) $(CFLAGS) -c ../av_mag.c

test_magvar.o : test_magvar.c ../av_mag.h

check : test_magvar
	./test_magvar WMM2020_TEST_VALUES.txt WMM2025_TestValues.txt
	python3 test_values.py WMM2020_TEST_VALUES.txt WMM2025_TestValues.txt
 
clean :
	- rm test_magvar test_magvar.o av_mag.o wmm.o

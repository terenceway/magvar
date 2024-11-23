# Build wmm.c (static coefficient data) given NGDC-provided COF
# files
#
# Example:
#	make -f wmm.mk

.PHONY : clean

wmm.c : compile-cof WMM2010.COF WMM2015.COF WMM2020.COF
	./compile-cof WMM2010.COF WMM2015.COF WMM2020.COF >wmm.c

clean :
	- rm wmm.c

include ../make.inc

all:  testdlamch testslamch testtimer install.out

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(LOADER) -o testdlamch dlamch.o lsame.o dlamchtst.o

testslamch: slamch.o lsame.o slamchtst.o
	$(LOADER) -o testslamch slamch.o lsame.o slamchtst.o

testtimer: superlu_timer.o timertst.o
	$(LOADER) -o testtimer superlu_timer.o timertst.o

install.out: install.csh
	@echo Testing machines parameters and timer 
	csh install.csh

slamch.o: slamch.c ; $(CC) -c $(NOOPTS) $<
dlamch.o: dlamch.c ; $(CC) -c $(NOOPTS) $<
superlu_timer.o: superlu_timer.c; $(CC) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o test* *.out

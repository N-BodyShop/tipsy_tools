TOOLS = ascii2bin bin2ascii totipnat totipstd snapshot
DEFS =  -DHAVE_STRING_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_ALLOCA=1
X_CFLAGS =   -nostdinc -I/usr/include -DSYSV -DSVR4 -DFUNCPROTO=7 -DNARROWPROTO
CFLAGS =  $(DEFS) $(X_CFLAGS)
all: $(TOOLS)

ascii2bin: ascii2bin.o
	$(CC) $(CFLAGS) -o ascii2bin ascii2bin.o $(LIBS)

bin2ascii: bin2ascii.o
	$(CC) $(CFLAGS) -o bin2ascii bin2ascii.o $(LIBS)

totipnat: totipnat.o
	$(CC) $(CFLAGS) -o totipnat totipnat.o $(LIBS)

totipstd: totipstd.o
	$(CC) $(CFLAGS) -o totipstd totipstd.o $(LIBS)

snapshot: snapshot.o
	$(CC) $(CFLAGS) -o snapshot snapshot.o $(LIBS)

clean:
	rm -f *.o
	rm -f $(TOOLS)

tipsy_tools:
	cd ..; tar cvf - tools/*.c tools/*.h tools/Makefile \
	| compress > tipsy_tools.tar.Z

ascii2bin.o: /usr/include/stdio.h tipsydefs.h
bin2ascii.o: /usr/include/stdio.h tipsydefs.h
totipnat.o: /usr/include/stdio.h tipsydefs.h /usr/include/rpc/types.h
totipnat.o: /usr/include/sys/types.h /usr/include/sys/time.h
totipnat.o: /usr/include/sys/time.h /usr/include/stdlib.h
totipnat.o: /usr/include/rpc/xdr.h
totipstd.o: /usr/include/stdio.h tipsydefs.h /usr/include/rpc/types.h
totipstd.o: /usr/include/sys/types.h /usr/include/sys/time.h
totipstd.o: /usr/include/sys/time.h /usr/include/stdlib.h
totipstd.o: /usr/include/rpc/xdr.h

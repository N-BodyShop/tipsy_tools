TOOLS = ascii2bin bin2ascii totipnat totipstd snapshot bin2genga
CFLAGS =  -Wall -O -g -I /usr/include/tirpc
LIBS = -ltirpc
all: $(TOOLS)

ascii2bin: ascii2bin.o
	$(CC) $(CFLAGS) -o ascii2bin ascii2bin.o $(LIBS)

ascii2std: ascii2std.o
	$(CC) $(CFLAGS) -o ascii2std ascii2std.o $(LIBS)

bin2ascii: bin2ascii.o
	$(CC) $(CFLAGS) -o bin2ascii bin2ascii.o $(LIBS)

totipnat: totipnat.o
	$(CC) $(CFLAGS) -o totipnat totipnat.o $(LIBS)

totipstd: totipstd.o
	$(CC) $(CFLAGS) -o totipstd totipstd.o $(LIBS)

snapshot: snapshot.o
	$(CC) $(CFLAGS) -o snapshot snapshot.o $(LIBS)

tipsy2snap: tipsy2snap.o
	$(CC) $(CFLAGS) -o tipsy2snap tipsy2snap.o -lm $(LIBS)

bin2genga: bin2genga.o
	$(CC) $(CFLAGS) -o bin2genga bin2genga.o -lm $(LIBS)

snap2tipsy: snap2tipsy.o
	$(CC) $(CFLAGS) -o snap2tipsy snap2tipsy.o -lm $(LIBS)

trimstd: trimstd.o
	$(CC) $(CFLAGS) -o trimstd trimstd.o -lm $(LIBS)

clean:
	rm -f *.o
	rm -f $(TOOLS)

tipsy_tools:
	cd ..; tar cvf - tools/*.c tools/*.h tools/Makefile \
	| compress > tipsy_tools.tar.Z

ascii2bin.o: /usr/include/stdio.h tipsydefs.h
bin2ascii.o: /usr/include/stdio.h tipsydefs.h
totipnat.o: /usr/include/stdio.h tipsydefs.h
totipnat.o: /usr/include/sys/types.h /usr/include/sys/time.h
totipnat.o: /usr/include/sys/time.h /usr/include/stdlib.h
totipstd.o: /usr/include/stdio.h tipsydefs.h
totipstd.o: /usr/include/sys/types.h /usr/include/sys/time.h
totipstd.o: /usr/include/sys/time.h /usr/include/stdlib.h

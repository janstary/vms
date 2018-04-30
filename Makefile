PREFIX	= $(HOME)
BINDIR	= $(PREFIX)/bin
MANDIR	= $(PREFIX)/man/man1

YFLAGS	= -d
CFLAGS	= -std=c99 -Wall -pedantic -pipe -D_BSD_SOURCE
# Yes, the -D_BSD_SOURCE is needed for the broken GNU

PROG	= vms
BINS	= ort ori
MAN1	= vms.1 ort.1 ori.1

all: $(BINS) $(MAN1)

atoms.o: atoms.c atoms.h
	$(CC) $(CFLAGS) -c atoms.c

request.o: request.c request.h
	$(CC) $(CFLAGS) -c request.c

ort: ort.c request.o atoms.o parse.y
	$(YACC) $(YFLAGS) parse.y
	$(CC) $(CFLAGS) -o ort ort.c request.o atoms.o y.tab.c

ori: ori.c request.o parse.y
	$(YACC) $(YFLAGS) parse.y
	$(CC) $(CFLAGS) -o ori ori.c request.o y.tab.c

lint: $(MAN1)
	mandoc -Tlint $(MAN1)

test: $(BINS) $(PROG)
	./vms -j 2 db.xyz request.single   ./single
	./vms -j 2 db.xyz request.optimize ./optimize
	./vms -j 2 db.xyz request.excited  ./excited

install: $(BINS) $(MAN1) #test
	install -d $(BINDIR)
	install -d $(MANDIR)
	install $(PROG) $(BINDIR)/
	install $(BINS) $(BINDIR)/
	install $(MAN1) $(MANDIR)/

uninstall:
	cd $(BINDIR) && rm -f $(BINS) $(PROG)
	cd $(MANDIR) && rm -f $(MAN1)

clean:
	rm -f $(BINS) *.o parse.c *.core y.* lex.* *~
	rm -rf ./single ./optimize ./excited

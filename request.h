#ifndef _REQUEST_H_
#define _REQUEST_H_

#include <stdio.h>

struct multi {
	unsigned even;
	unsigned odd;
};

struct absorb {
	float  lolen;
	float  hilen;
	float  loint;
	float  hiint;	/* < 0 means unlimited */
	float  best;	/* -1 = max, -2 = mid */
	float  weight;
	double maxeval;
};

struct request {
	enum {	NOJOB,
		SINGLE,
		OPTIMIZE,
		EXCITED
	}		job;
	char*		basis;
	char*		method;
	int		charge;
	struct multi*	multi;
	unsigned	memory;
	unsigned	nroots;
	unsigned	iroot;
	unsigned	triplets;
	unsigned	maxdim;
	unsigned	nbest;

	struct absorb**	absorb;
	unsigned	numabsorb;
};

struct request* parsereq(FILE*);
int reqcheck(struct request*);
void printreq(struct request*);
void freereq(struct request*);

#endif

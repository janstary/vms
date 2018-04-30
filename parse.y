%{
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <err.h>

#include "request.h"

struct request *request;

static void yyerror(char *);
static int yylex(void);
FILE *yyfp;

typedef struct {
	union {
		int   i;
		float f;
		char* s;
	};
	int lineno;
	int col;
} yystype;
#define YYSTYPE yystype

%}
%token JOB JOBNAME BASIS METHOD MEMORY
%token CHARGE MULTI TRIPLETS NROOTS IROOT MAXDIM NBEST ABSORB
%token STRING INTEGER UNSIGNED FLOAT

%%
input	: /* empty */
	| input		'\n'
	| input line	'\n'
	| input error	'\n'
		{ yyerrok; /* FIXME not OK */ }
	;

line	: JOB    JOBNAME		{ request->job    = $2.i; }
	| BASIS  STRING 		{ request->basis  = $2.s; }
	| METHOD STRING			{ request->method = $2.s; }
	| MEMORY UNSIGNED		{ request->memory = $2.f; }
	| CHARGE integer		{ request->charge = $2.f; }
	| MULTI  UNSIGNED UNSIGNED	{
		if (request->multi) {
			warnx("More than one multiplicity directive");
			/* FIXME we should probably error out here */
			free(request->multi);
		}
		if ((request->multi = calloc(1, sizeof(struct multi))) == NULL)
			err(1, NULL);
		request->multi->even = $2.f;
		request->multi->odd  = $3.f;
	}
	| TRIPLETS			{ request->triplets = 1; }
	| NROOTS UNSIGNED		{ request->nroots = $2.f; }
	| IROOT UNSIGNED		{ request->iroot  = $2.f; }
	| MAXDIM UNSIGNED		{ request->maxdim = $2.f; }
	| NBEST  UNSIGNED		{ request->nbest  = $2.f; }
	| ABSORB number number number number best {
		struct absorb *a;
		if ((a = calloc(1, sizeof(*a))) == NULL)
			err(1, "Cannot allocate absorb");
		a->lolen = $2.f;
		a->hilen = $3.f;
		a->loint = $4.f;
		a->hiint = $5.f;
		a->best  = $6.f;
		if ((request->absorb = realloc(request->absorb,
		sizeof(a) * (1 + request->numabsorb))) == NULL)
			err(1, "Cannot allocate absorbs");
		request->absorb[request->numabsorb++] = a;
	}
	;

integer	: INTEGER
	| UNSIGNED
	;

number	: integer
	| FLOAT
	;

best	: { $$.f = -2; }
	| number
	;
%%

void
yyerror(char *msg) {
	warnx("%s", msg);
}

static struct {
	const char *word;
	int token;
} keywords[] = {
	{ "job",	JOB    },
	{ "basis",	BASIS  },
	{ "method",	METHOD },
	{ "memory",	MEMORY },
	{ "charge",	CHARGE },
	{ "multi",	MULTI },
	{ "triplets",	TRIPLETS },
	{ "nroots",	NROOTS },
	{ "iroot",	IROOT },
	{ "maxdim",	MAXDIM },
	{ "nbest",	NBEST  },
	{ "absorb",	ABSORB },
};

static struct {
	const char *word;
	int job;
} jobs[] = {
	{ "single",	SINGLE },
	{ "optimize",	OPTIMIZE },
	{ "excited",	EXCITED },
};

int
yylex(void)
{
	int i, c;
	char *str, *e;
	char buf[1024], *p = buf;
	char *ebuf = buf + sizeof(buf);
	float f;

repeat:
	while ((c = getc(yyfp)) == ' ' || c == '\t')
		;

	if (c == '\n')
		return c;
	if (c == EOF)
		return 0;

	for (;; c = getc(yyfp)) {
		switch (c) {
			case EOF:
				return 0;
				break;
			case '#' :
				while ((c = getc(yyfp)) != '\n')
					;
				/* FALLTRHOUGH */
			case ' ' :
			case '\t':
			case '\n':
				goto eow;
				break;
		}
		*p++ = c;
		if (p == ebuf) {
			yyerror("line too long");
			p = buf;
		}
	}

eow:
	*p = 0;
	if (c != EOF)
		ungetc(c, yyfp);
	if (p == buf) {
		goto repeat;
	}

	for (i = 0; i < sizeof(keywords) / sizeof(keywords[0]); i++) {
		if (0 == strcasecmp(buf, keywords[i].word))
			return keywords[i].token;
	}

	for (i = 0; i < sizeof(jobs) / sizeof(jobs[0]); i++) {
		if (0 == strcasecmp(buf, jobs[i].word)) {
			yylval.i = jobs[i].job;
			return JOBNAME;
		}
	}

	i = strtol(buf, &e, 10);
	if (*e == 0) {
		yylval.f = i;
		return i > 0 ? UNSIGNED : INTEGER;
	}

	f = strtof(buf, &e);
	if (*e == 0) {
		yylval.f = f;
		return FLOAT;
	}

	if ((str = strdup(buf)) == NULL)
		err(1, "strdup");
	yylval.s = str;
	return STRING;
}

struct request*
parsereq(FILE *in)
{
	yyfp = in;
	request = calloc(1, sizeof(struct request));
	return (request && (yyparse() == 0)) ? request : NULL;
}

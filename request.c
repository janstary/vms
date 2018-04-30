#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <err.h>

#include "request.h"

int
reqcheck(struct request* r)
{
	int i;

	if (r->method == NULL) {
		warnx("Method not set.");
		return -1;
	}
	if (r->memory == 0) {
		r->memory = 1024;
	}
	
	switch (r->job) {
	case NOJOB:
		warnx("Job not specified in request.");
		return -1;
		break;
	case SINGLE:
		break;
	case OPTIMIZE:
		if (0 == r->maxdim) {
			r->maxdim = 10 * r->iroot;
		}
		break;
	case EXCITED:
		if (0 == r->nroots) {
			warnx("nroots not set for an excited job");
			return -1;
		}
		if (0 == r->maxdim) {
			r->maxdim = 10 * r->nroots;
		}
		if (r->multi) {
			if ((r->multi->even % 2 == 0)
			||  (r->multi->odd  % 2 != 0)) {
				warnx("wrong multiplicity %u %u",
					r->multi->even, r->multi->odd);
				return -1;
			}
		}
		for (i = 0; i < r->numabsorb; i++) {
			struct absorb *a = r->absorb[i];
			if (a->lolen < 0 || a->hilen < 0) {
				warnx("Wavelength %f to %f?",a->lolen,a->hilen);
				return -1;
			}
			if (a->lolen > a->hilen) {
				warnx("lo > hi: %f > %f", a->lolen, a->hilen);
				return -1;
			}
			if (a->loint < 0) {
				warnx("Low intesity negative: %f", a->loint);
				return -1;
			}
			if (a->hiint > 0
			&&  a->loint > a->hiint) {
				warnx("lo > hi: %f > %f", a->loint, a->hiint);
				return -1;
			}
			if (!(a->best == -1    /* max is best */
			  ||  a->best == -2    /* mid is best */
			  ||  a->best >= 0)) { /* the given number is best */
				warnx("best value %f makes no sense", a->best);
				return -1;
			}
		}
		break;
	}
	return 0;
}

void
printreq(struct request* r)
{
	int i;
	struct absorb *a;

	if (r == NULL || r->job == 0)
		return;

	printf("# User request translated for ORCA by ort(1).\n");
	printf("# Do not edit this, edit the original request.\n#\n");

	switch (r->job)
	{
		case SINGLE:
			printf("# single point calculations miniprint\n");
			printf("\n!%s %s\n\n", r->method,
				r->basis ? r->basis : "");
			printf("%%maxcore %u\n", r->memory);
			break;
		case OPTIMIZE:
			printf("# geometry optimization\n");
			printf("\n!%s %s Opt miniprint\n", r->method,
				r->basis ? r->basis : "");
			printf("%%maxcore %u\n", r->memory);
			if (r->iroot) {
				printf(
					"%%tddft\n"
					" nroots %u\n"
					" iroot  %u\n"
					" maxdim %u\n"
					" end\n\n",
					r->iroot, r->iroot, r->maxdim);
			}
			break;
		case EXCITED:
			printf("# computing excited states\n");
			if (r->charge) {
				printf("# charge %d", r->charge);
				if (r->multi)
					printf(", multiplicity %u %u",
						r->multi->even, r->multi->odd);
				putchar('\n');
			}
			for (i = 0; i < r->numabsorb && (a=r->absorb[i]); i++) {
				printf("# absorb %09.5f %09.5f %09.5f %09.5f ",
					a->lolen, a->hilen, a->loint, a->hiint);
				if (a->best == -2) {
					printf("[middle]\n");
				} else if (a->best == -1) {
					printf("[maximal]\n");
				} else {
					printf("%f\n", a->best);
				}
			}
			printf(
				"\n!%s %s miniprint\n"
				"%%maxcore %u\n"
				"%%tddft\n"
				" nroots %u\n"
				" maxdim %u\n"
				" %s"
				" end\n\n",
				r->method, r->basis ? r->basis : "",
				r->memory, r->nroots, r->maxdim,
				r->triplets ? "triplets true\n" : "");
			break;
		default:
			break;
	}
}

void
freereq(struct request* r)
{
	int i;
	if (r) {
		if (r->absorb) {
			for (i = 0; r->numabsorb; r->numabsorb--, i++)
				free(r->absorb[i]);
			free(r->absorb);
		}
		free(r);
	}
}

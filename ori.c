#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <fts.h>
#include <err.h>

#include "request.h"


extern char *optarg;
extern int opterr;
extern int optind;

int g_quiet = 0;

#define WARN(...) while(g_quiet == 0 && (warn(__VA_ARGS__), 0))
#define WARNX(...) while(g_quiet == 0 && (warnx(__VA_ARGS__), 0))

#define REQNAME "request"
#define MOLNAME "molecule."
#define GOODNAME(name) (strstr(name, MOLNAME) == name)

struct atom {
	char	*name;		/* "He" etc. */
	float	x, y, z;	/* position */
};

struct state {
	unsigned	num;	/* order */
	float		wlen;	/* wavelength */
	float		inty;	/* intensity of transition */
};

struct molecule {
	unsigned	id;
	char*		filename;
	unsigned	selected;
	unsigned	numatoms;
	struct atom	**atoms;
	unsigned	numstates;
	struct state	**states;
	double		eval;
};

struct request *request = NULL;
struct molecule **molecules;
unsigned sizemols = 1;
unsigned havemols = 0;

void
usage()
{
	extern char* __progname;
	fprintf(stderr, "usage: %s orcadir\n", __progname);
}

int
is_blank_line(const char* line)
{
	const char *p;
	for (p = line; isspace(*p); p++)
		;
	return *p == '\0';
}

/* Look for a line in an open file, advancing the pointer.
 * Skip ' ' and also '*' as those are mostly ORCA decorations.
 * Return 1 if found, 0 if not. */
int
grep(const char* hdr, FILE* file)
{
	char *line = NULL;
	size_t size = 0;
	ssize_t len;
	char *p;
	while ((len = getline(&line, &size, file)) != -1) {
		for (p = line; *p == ' ' || *p == '*' ; p++)
			;
		if (0 == strncmp(p, hdr, strlen(hdr)))
			return 1;
	}
	return 0;
}

/* Assign id to molecule from its filename */
unsigned
molnum(char* filename)
{
	unsigned id = 0;
	char *p, *err = NULL;
	if ((filename == NULL) || !GOODNAME(filename))
		return 0;
	p = strrchr(filename, '.');
	id = strtol(++p, &err, 10);
	return *err ? 0 : id;
}

struct molecule*
mkmol(FTSENT *molfile)
{
	struct molecule *mol;
	if (molfile == NULL)
		return NULL;
	if ((mol = calloc(1, sizeof(struct molecule))) == NULL)
		err(1, NULL);
	mol->filename = strdup(molfile->fts_path);
	if (0 == (mol->id = molnum(molfile->fts_name)))
		WARNX("No number in '%s'", mol->filename);
	return mol;
}

unsigned
addmol(struct molecule *mol) {
	if (mol == NULL)
		return -1;
	if (havemols == sizemols) {
		sizemols *= 2;
		molecules=realloc(molecules,sizemols*sizeof(struct molecule*));
		if (molecules == NULL)
			err(1, "Cannot reallocate for %u molecules", sizemols);
	}
	molecules[havemols++] = mol;
	return -1;
}

void
pratom(struct atom *a)
{
	if (a)
		printf("%s %9.6f %9.6f %9.6f\n", a->name, a->x, a->y, a->z);
}

void
prstate(struct state *s)
{
	if (s)
		fprintf(stderr, "%03d %f %f\n", s->num, s->wlen, s->inty);
}

void
prmol_single(struct molecule *mol)
{
	printf("molecule %u from %s\n", mol->id, mol->filename);
}

void
prmol_optimize(struct molecule *mol)
{
	unsigned i;
	if (mol == NULL)
		return;
	printf("%u\n%s\n", mol->numatoms, mol->filename);
	for (i = 0; i < mol->numatoms; i++)
		pratom(mol->atoms[i]);
}

void
prmol_excited(struct molecule *mol)
{
	unsigned i;
	if (mol == NULL)
		return;
	printf("%u\n%s\n", mol->numatoms, mol->filename);
	for (i = 0; i < mol->numatoms; i++)
		pratom(mol->atoms[i]);
	if (g_quiet == 1) {
		/* Be concise so that mm(1) and others can parse our output. */
		fprintf(stderr, "%u\t%f\n", mol->id, mol->eval);
	} else {
		fprintf(stderr, "molecule %u from %s", mol->id, mol->filename);
		if (mol->eval)
		fprintf(stderr, ": eval %f", mol->eval);
		fprintf(stderr, "\n");
		for (i = 0; i < mol->numstates; i++)
			prstate(mol->states[i]);
	}
}

void
print_single()
{
	unsigned i;
	if (molecules == NULL)
		return;
	for (i = 0; i < havemols; i++)
		prmol_single(molecules[i]);
}
 
void
print_optimize()
{
	unsigned i;
	if (molecules == NULL)
		return;
	for (i = 0; i < havemols; i++)
		prmol_optimize(molecules[i]);
}

void
print_excited(int selected)
{
	unsigned i, max;
	if (molecules == NULL)
		return;
	if (request->nbest == 0 || request->nbest > havemols)
		max = havemols;
	else
		max = request->nbest;
	for (i = 0; i < max; i++)
		if (selected == 0 || molecules[i]->selected)
			prmol_excited(molecules[i]);
}

/* Single point calculations: nothing to evaluate. */
int
ori_single(FTS *dir)
{
	return 0;
}

/* Read one line of geometry, such as
   O  0.101077   0.000000    0.000000
   H  0.708064   0.000000    0.771049
   H  0.708064   0.000000   -0.771049
*/
int
read_optimize_line(struct molecule *mol, char* line)
{
	char *p;
	unsigned n;
	struct atom *a;
	if (mol == NULL)
		return -1;
	if ((a = calloc(1, sizeof(struct atom))) == NULL)
		err(1, NULL);
	n = mol->numatoms;
	if ((mol->atoms = realloc(mol->atoms,
	(1+n) * sizeof(struct atom*))) == NULL)
		err(1, NULL);
	do { p = strsep(&line, " "); } while (*p == '\0');
	a->name = strdup(p);
	do { p = strsep(&line, " "); } while (*p == '\0');
	a->x = strtof(p, NULL);
	do { p = strsep(&line, " "); } while (*p == '\0');
	a->y = strtof(p, NULL);
	do { p = strsep(&line, " "); } while (*p == '\0');
	a->z = strtof(p, NULL);
	mol->atoms[mol->numatoms++] = a;
	return 0;
}

/* Read the optimized coordinate lines from an ORCA result file
 * and save the optimized coordinates in the molecule structure. */
int
read_optimize(struct molecule *mol, char* orcalog)
{
	FILE *orca;
	char *line = NULL;
	size_t size = 0;
	ssize_t len;
	const char *hdr;
	if ((orca = fopen(orcalog, "r")) == NULL) {
		WARN("Cannot open '%s'", orcalog);
		return -1;
	}

	hdr = "THE OPTIMIZATION HAS CONVERGED";
	if (!grep(hdr, orca)) {
		WARNX("'%s'\nnot found in %s", hdr, orcalog);
		goto bad;
	}
	hdr = "CARTESIAN COORDINATES (ANGSTROEM)";
	if (!grep(hdr, orca)) {
		WARNX("'%s'\nnot found in %s", hdr, orcalog);
		goto bad;
	}
	if (!grep("---------", orca)) {
		WARNX("delimiter lines not found");
		goto bad;
	}
	while ((len = getline(&line, &size, orca)) != -1) {
		if (is_blank_line(line))
			break;
		read_optimize_line(mol, line);
	}
	free(line);
	fclose(orca);
	return mol->numatoms;
bad:
	free(line);
	fclose(orca);
	return 0;
}

/* Geometry optimizations: nothing to evaluate,
 * just print out the resulting geometry. */
int
ori_optimize(FTS *dir)
{
	FTSENT *molfile;
	struct molecule *mol;
	if ((molecules = calloc(sizemols, sizeof(struct molecule*))) == NULL)
		err(1, NULL);
	while ((molfile = fts_read(dir))) {
		if ((molfile->fts_info != FTS_F)
		|| !GOODNAME(molfile->fts_name))
			continue;
		mol = mkmol(molfile);
		if (!read_optimize(mol, molfile->fts_path))
			WARNX("no atoms in '%s'",
				mol->filename);
		if (mol->numstates < request->nroots)
			WARNX("%u < %u states in '%s'",
				mol->numstates, request->nroots, mol->filename);
		addmol(mol);
	}
	print_optimize();
	return 0;
}

/* Read one excited state, such as
State   Energy  Wavelength   fosc         T2         TX        TY        TZ
   1   10976.2    911.1   0.012410144   0.37222   0.60429  -0.07878  -0.02908
*/
int
read_excited_line(struct molecule *mol, char* line)
{
	char *p;
	unsigned iroot;
	struct state *s;
	if (mol == NULL)
		return -1;
	do { p = strsep(&line, " "); } while (*p == '\0');
	iroot = strtoul(p, NULL, 10);
	if (request->iroot && iroot != request->iroot){
		return 0;
	}
	if ((s = calloc(1, sizeof(struct state))) == NULL)
		err(1, NULL);
	if ((mol->states = realloc(mol->states,
	(mol->numstates + 1) * sizeof(struct state*))) == NULL)
		err(1, NULL);
	s->num = iroot;
	do { p = strsep(&line, " "); } while (*p == '\0');
	do { p = strsep(&line, " "); } while (*p == '\0');
	s->wlen = strtof(p, NULL);
	do { p = strsep(&line, " "); } while (*p == '\0');
	if (strcmp(p,"spin") == 0)
		s->inty = -1;
	else
		s->inty = strtof(p, NULL);
	mol->states[mol->numstates++] = s;
	return 0;
}

/* Read the excited state lines from an ORCA result file
 * and save the excited states in the molecule structure.
 * FIXME Return the number of excited states, or -1 on error.*/
int
read_excited(struct molecule *mol, char* orcalog)
{
	FILE *orca;
	char *line = NULL;
	size_t size = 0;
	ssize_t len;
	const char *hdr;
	if ((orca = fopen(orcalog, "r")) == NULL) {
		WARN("Cannot open '%s'", orcalog);
		return -1;
	}
	hdr = "CARTESIAN COORDINATES (ANGSTROEM)";
	if (!grep(hdr, orca)) {
		WARNX("'%s'\nnot found in %s", hdr, orcalog);
		goto bad;
	}
	if (!grep("---------", orca)) {
		WARNX("delimiter lines not found");
		goto bad;
	}
	while ((len = getline(&line, &size, orca)) != -1) {
		if (is_blank_line(line))
			break;
		read_optimize_line(mol, line);
	}
	hdr = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS";
	if (!grep(hdr, orca)) {
		WARNX("'%s'\nnot found in %s", hdr, orcalog);
		goto bad;
	}
	if (!grep("--------", orca)) {
		WARNX("delimiter not found in %s", orcalog);
		goto bad;
	}
	if (!grep("--------", orca)) {
		WARNX("delimiter not found in %s", orcalog);
		goto bad;
	}
	while ((len = getline(&line, &size, orca)) != -1) {
		if (is_blank_line(line))
			break;
		read_excited_line(mol, line);
	}
	free(line);
	fclose(orca);
	return mol->numstates;
bad:
	free(line);
	fclose(orca);
	return 0;
}

/* Discover the molecules which satisfy the absorption criteria.
 * Return the number of molecules passing the inequalities. */
unsigned
filter_absorb()
{
	struct molecule *mol;
	struct absorb *abs;
	struct state *sta;
	int m, a, s, good = 0;

	if (request->absorb == NULL)
		WARNX("No absorption criteria. All molecules pass.");

	for (m = 0; (m < havemols) && (mol = molecules[m]); m++) {
		for (a = 0; a < request->numabsorb; a++) {
			abs = request->absorb[a];
			mol->selected = 0;
			for (s = 0; s < mol->numstates; s++) {
				sta = mol->states[s];
				if (sta->wlen >= abs->lolen
				&&  sta->wlen <= abs->hilen
				&&  sta->inty >= abs->loint
				&& (sta->inty <= abs->hiint || abs->hiint < 0)){
					mol->selected = 1;
					break;
				}
			}
			if (mol->selected == 0)
				break;
		}
		good += mol->selected;
        }
	return good;
}

/* Calculate the maximal intensity
 * specified in an absorb filter */
float
maxint(struct absorb *abs)
{
	struct molecule *mol;
	struct state *sta;
	int m, s = 0;
	float maxint = abs->loint;
	for (m = 0; (m < havemols) && (mol = molecules[m]); m++) {
		if (mol->selected != 1)
			continue;
		for (s = 0; s < mol->numstates; s++) {
			sta = mol->states[s];
			if (sta->wlen >= abs->lolen
			&&  sta->wlen <= abs->hilen
			&&  sta->inty > maxint) {
				maxint = sta->inty;
			}
		}
	}
	return maxint;
}

int compare(const void * a, const void * b)
{
	struct molecule *molA = *(struct molecule **)a;
	struct molecule *molB = *(struct molecule **)b;
	return
		(molA->eval == molB->eval) ? 0 :
		(molA->eval < molB->eval)  ? 1 : -1;
}

void
sort_absorb()
{

	struct molecule *mol;
	struct absorb *abs;
	struct state *sta;
	int m, a, s = 0;
	double eval, besteval;

	if (request->absorb == NULL)
		WARNX("No absorption criteria. All molecules pass.");

	for (a = 0; a < request->numabsorb; a++) {
		abs = request->absorb[a];
		if (abs->hiint == -1)
			abs->hiint = maxint(abs);
		if (abs->best == -1) /* max */
			abs->best = abs->hiint;
		else if (abs->best == -2) /* mid */
			abs->best = (abs->hiint + abs->loint) / 2;
		if(abs->hilen != abs->lolen
		&& abs->hiint != abs->loint)
			abs->weight = 1 / (abs->hilen - abs->lolen)
					/ (abs->hiint - abs->loint);
		else
			abs->weight = 1000000;
		abs->maxeval =
			(abs->best - abs->loint < abs->hiint - abs->best)
			? abs->hiint - abs->best
			: abs->best - abs->loint;
		if(!abs->maxeval)
			abs->maxeval = 1;
        }

	for (m = 0; (m < havemols) && (mol = molecules[m]); m++) {
		if (mol->selected != 1)
			continue;
		for (a = 0; a < request->numabsorb; a++) {
			abs = request->absorb[a];
			besteval = abs->maxeval;
			for (s = 0; s < mol->numstates; s++) {
				sta = mol->states[s];
				if (sta->wlen >= abs->lolen
				&&  sta->wlen <= abs->hilen
				&&  sta->inty >= abs->loint
				&&  sta->inty <= abs->hiint) {
					eval = fabs(abs->best - sta->inty);
					if (eval < besteval){
						besteval = eval;
					}
				}
			}
			mol->eval += (1 - besteval/abs->maxeval) * abs->weight;
		}
        }
	qsort(molecules, havemols, sizeof(struct molecule*), compare);
}

/* Excited state calculations:
 * evaluate according to request. */
int
ori_excited(FTS *dir)
{
	FTSENT *molfile;
	struct molecule *mol;
	if ((molecules = calloc(sizemols, sizeof(struct molecule*))) == NULL)
		err(1, NULL);
	while ((molfile = fts_read(dir))) {
		if ((molfile->fts_info != FTS_F)
		|| !GOODNAME(molfile->fts_name))
			continue;
		mol = mkmol(molfile);
		if (!read_excited(mol, molfile->fts_path))
			WARNX("no excited states in '%s'",
				mol->filename);
		if (mol->numstates < request->nroots && !request->iroot)
			WARNX("number of states %u < %u states in '%s'",
				mol->numstates, request->nroots, mol->filename);
		else if (request->iroot && mol->numstates != 1)
			WARNX("number of states %u != 1 state for root %u in '%s'",
				mol->numstates, request->iroot, mol->filename);
		addmol(mol);
	}
	if (request->absorb) {
		filter_absorb();
		if (request->nbest || g_quiet)
			sort_absorb();
		print_excited(1);
	} else {
		print_excited(0);
	}
	return 0;
}

int
main(int argc, char** argv)
{
	char  reqname[1024];
	FILE* reqfile;


	FTS *dir = NULL;
	char* const *orcadir;

	int c;

	while ((c = getopt(argc, argv, "q")) != -1) switch (c) {
	case 'q':
		g_quiet = 1;
		break;
	default:
		usage();
		return 1;
	}
	argc -= optind;
	argv += optind;

	if (argc != 1) {
		usage();
		return 1;
	}

	snprintf(reqname, 1024, "%s/%s", argv[0], REQNAME);
	if ((reqfile = fopen(reqname, "r")) == NULL) {
                WARN("Cannot fopen %s", reqname);
                return 1;
	}
        if ((request = parsereq(reqfile)) == NULL) {
                WARNX("Cannot parse the request file.");
                return 1;
        }
        if (reqcheck(request) != 0) {
                WARNX("The request is not valid.");
                return 1;
        }

	orcadir = &argv[0];
	if ((dir = fts_open(orcadir, FTS_LOGICAL|FTS_NOCHDIR|FTS_XDEV, NULL))
	== NULL)
		err(1, "Cannot traverse %s", *orcadir);
	switch (request->job) {
		case NOJOB:
			break;
		case SINGLE:
			ori_single(dir);
			break;
		case OPTIMIZE:
			ori_optimize(dir);
			break;
		case EXCITED:
			ori_excited(dir);
			break;
		default:
			break;
	}
	fts_close(dir);
        freereq(request);
	return 0;
}

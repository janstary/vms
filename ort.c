#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <err.h>

#include "request.h"
#include "atoms.h"

void
usage()
{
	extern char* __progname;
	fprintf(stderr, "usage: %s request [molecule]\n", __progname);
}

/* Figure out the charge-an-multiplicity line,
 * depending on the molecule at hand and the request. */
int
cm(FILE* file, unsigned charge, struct multi *multi)
{
	size_t size;
	ssize_t len;
	char *line = NULL;
	char *atom = NULL;

	unsigned int ele = 0;
	unsigned int sum = 0;
	unsigned int multiplicity = 0;

	while ((len = getline(&line, &size, file)) != -1) {
		atom = strsep(&line, " \t");
		if (line == NULL)
			return -1;
		line = atom;
		if ((ele = getatom(atom)) == -1)
			return -1;
		sum += ele;
	}
	free(line);
	if (charge)
		sum -= charge;
	multiplicity = (sum % 2)
		? (multi && multi->odd  > 2 ? multi->odd : 2)
		: (multi && multi->even > 1 ? multi->even : 1);
	printf("%d %d\n", charge, multiplicity);
	return 0;
}

int
main(int argc, char** argv)
{
	FILE* reqfile;
	FILE* molfile;
	struct request *request;

	if (argc < 2 || argc > 3) {
		usage();
		return 1;
	}

	if ((reqfile = fopen(argv[1], "r")) == NULL)
		err(1, "Cannot fopen '%s'", argv[1]);
	if ((request = parsereq(reqfile)) == NULL) {
		warnx("Cannot parse '%s'", argv[1]);
		return 1;
	}
	if (reqcheck(request) != 0) {
		warnx("Request in '%s' is not valid", argv[1]);
		return 1;
	}

	if (argc == 2)
		printreq(request);
	if (argc == 3) {
		if ((molfile = fopen(argv[2], "r")) == NULL)
			err(1, "Cannot fopen '%s'", argv[2]);
		cm(molfile, request->charge, request->multi);
		fclose(molfile);
	}

	fclose(reqfile);
	freereq(request);
	return 0;
}

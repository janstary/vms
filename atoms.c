#include <strings.h>

#define La "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"
#define Ac "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"

/* The index into this array
 * is the actual atomic number. */
const char* atoms[] = {
	"Adamant",
	"H",						  							 "He",
	"Li", "Be",								"B",   "C",  "N",   "O",  "F",   "Ne",
	"Na", "Mg",								"Al",  "Si", "P",   "S",  "Cl",  "Ar",
	"K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",  "Ge", "As",  "Se", "Br",  "Kr",
	"Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",  "Sn", "Sb",  "Te", "I",   "Xe",
	"Cs", "Ba",  La,  "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",  "Pb", "Bi",  "Po", "At",  "Rn",
	"Fr", "Ra",  Ac,  "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl", "Uup", "Lv", "Uus", "Uuo",
};

#define NUMATOMS (sizeof(atoms) / sizeof(const char*))

int
getatom(const char* name)
{
	int i = -1;
	for (i = 0; i < NUMATOMS; i++)
		if (0 == strcasecmp(name, atoms[i]))
			return i;
	return -1;
}

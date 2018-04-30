# Virtual Molecule Screener

VMS builds upon the
[ORCA](https://orcaforum.cec.mpg.de/) quantum chemistry software.
It parses a user's request (such as: light absorption),
submits the appropriate instructions to ORCA,
along with a database of molecules,
and interprets the results.

The application is to test moleules "in silico"
before actually going to the lab, spending precious time and money.

## installation

You will need to have [ORCA](https://orcaforum.cec.mpg.de/) installed
somewhere in your `PATH`, obviously.
Currently, `vms` only supports ORCA 3,
and needs to be updated to use ORCA 4.

To install vms itself,

```sh
mkdir -p ~/src && cd ~/src
git clone git@github.com:janstary/vms.git
cd vms && make clean install
```

This will install `vms` into your `$HOME`.
You can tweak that (and other stuff) in the `Makefile`.
VMS should compile on any POSIX-compliant UNIX system.

## example

VMS comes with `db.xyz` containing just two molecules:
the molecule of water, and the molecule of methan.
An example run of `vms` can be seen with `make test`:

```sh
./vms -j 2 db.xyz request.single   ./single
./vms -j 2 db.xyz request.optimize ./optimize
./vms -j 2 db.xyz request.excited  ./excited
```

The `single` request instructs ORCA to run the relatively cheap
single point calculations on the molecules. This is a test that
`vms` basicaly works. Look into the resulting `./single` directory.

The `optimize` request instructs ORCA to optimize the molecule geometry.
Note the modified geometry in `./optimize/output`, optimized from the
naive structure of H2O and CH4 in `db.xyz`.

The `excited` request instructs ORCA to compute excited states,
in order to figure out the absorption of light:

```
absorb 300 500 0.100 -1.00 -1
```

This means we are looking for molecules that absorb light
with wavelenghts between 300 and 500 nm, with intesity at least 0.10
(which means, roughly, that a photon has
at least 10% chance of being absorbed);
the higher the absorption the better.

In `./excited/report`, you will see that

```
molecule 2 from ./excited/molecule.000002: eval 0.039725
003 385.399994 0.127865
```

This means that the methan molecule, in its third excited state,
absorbs light with a wavelength of 385 nm (violet) with an intensity of 12%,
therefore passing the specification in the request. The molecule of water
does not pass the filter; so `./excited/output` only contains the methan.

## documentation

VMS installs three binaries: `vms`, `ort` and `ori`.
Their manpages describe the details of usage,
the format of the request files, and the structure of the output.

## authors

VMS was originaly written by
[Štěpán Sršeň](http://photox.vscht.cz/people.php?person=srsen)
and later modified and extended by
[Jonatan Matějka](https://github.com/jonatan1024)
and Jan Starý.

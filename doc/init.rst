Initialization
====================

The main script is `init_1cpn.py`

It takes some command line arguments:

* `--stemangle` (`-a`) : This is the angle alpha
* `--nrl` (`-nrl`) : Nucleosome Repeat Length
* `--nrlends` (`-nrlends`) : Nucleosome Repeat Length of DNA at beginning and end of fiber. This arguement defaults to set `nrlends=nrl` if not specified
* `--nnucl` (`-n`): number of nucleosomes
* `--linkerhistone` (`-lh`): turn on linker histones



note that if -n = 0, then only DNA will be generated of length -nrl


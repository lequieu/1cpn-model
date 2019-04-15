.. _label-initialization:

Initialization
====================

Initial configurations for 1CPN can be generated using `input/init_1cpn.py`

`init_1cpn.py` can take various command line arguements:

* ``--stemangle`` (``-a``) : This is the angle alpha
* ``--nrl`` (``-nrl``) : Nucleosome Repeat Length
* ``--nrlends`` (``-nrlends``) : Nucleosome Repeat Length of DNA at beginning and end of fiber. This arguement defaults to set ``nrlends=nrl`` if not specified
* ``--nnucl`` (``-n``): number of nucleosomes
* ``--linkerhistone`` (``-lh``): turn on linker histones

.. note::
    If you only want to generate a long strand of DNA (without any nucleosomes), you can do this by adding ``--nnucl 0`` to the command line arguements. Then specify the length of DNA (in base pairs) with the ``--nrl`` entry.

.. note::
    If you use the ``-lh`` flag to turn generate the linker histone, you'll also  need to change a single line in ``in.1cpn``. You'll want to  change this line to 
    :: 
        variable lh equal 1






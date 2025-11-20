## CP2K-VAB
This is a working version of https://github.com/dupuislab/CP2K. Starting from the same commit as the
former link for version 7.1 (commit # af259e795f), relevent file changes were applied. Some additional
modifications were necessary for the code to compile with a more modern toolchain and ARCH file.

# How to use

Git clone this branch:
`git clone --recursive --branch cp2k-vab https://github.com/abussy/cp2k.git`

Ensure DBCSR version is compatible:
`cd exts/dbcsr/ && git checkout v2.0.0`

Adapt the `arch/example.psmp` ARCH file to match your local toolchain installation. If you need Libxc,
then make sure to install a compatible version (current v7.0.0 does not work. v4.3.4 should)

Compile in the usual way: `make -j 4 ARCH=example VERSION=psmp`.

The VAB tests from https://github.com/dupuislab/CP2K/tree/master/MOA_CP2K seem to pass. Note that many
of them use erroneous absolute path for the potential and basis set files (which need fixing). It is
**highly** recommended to run at least some of these tests afer compilation.

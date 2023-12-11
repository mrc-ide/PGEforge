## Installing `DEploidIBD`

*Please note that the following installation instructions were derived from those described in the [github repository](https://github.com/mcveanlab/DEploid#installation).*

**TODO:**
- For now only includes installation instructions for MacOS

#### Step 1: Preparing for C++ compilation.
`DEploidIBD` is implemented in C++ and the code needs to be compiled to your local computer hardware before it can run.

**MacOS**
Use the package manager [Homebrew](https://brew.sh/) to install `autoconf` and other tools necessary for compilation:
```
brew install automake autoconf autoconf-archive cppunit
```

#### Step 2: Clone the repository.
Use [git](https://git-scm.com/) to clone the repository:
```
git clone https://github.com/DEploid-dev/DEploid.git
cd DEploid
git submodule update --init --recursive --remote
```

#### Step 3: Compile!
Finally, you can compile `DEploidIBD` with the following commands:
```
./bootstrap
make
```

If the compilation has run successfully, you should be able to run the executable with:
```
./dEploid
```

The github repository also includes some example data and commands. For example, you can run:

```
./dEploid \
 -vcf data/testData/PG0390-C.test.vcf \
 -plaf data/testData/labStrains.test.PLAF.txt 
 -o PG0390-CNopanel \
 -noPanel
```

which  will produce a set of example output files.


## Installing `DEploidIBD`

*Please note that the following installation instructions were derived from those described in the [github repository](https://github.com/mcveanlab/DEploid#installation).*

#### Step 1: Preparing for C++ compilation.
`DEploidIBD` is implemented in C++ and the code needs to be compiled to your local computer hardware before it can run.

**MacOS**
Use the package manager [Homebrew](https://brew.sh/) to install `autoconf` and other tools necessary for compilation:
```
brew install pkg-config automake autoconf autoconf-archive cppunit
```

**Ubuntu/Debian**
For Ubuntu or Debian operaing systems, use the package manager `apt-get` to install the compilation tools:
```
apt-get install build-essential autoconf autoconf-archive libcppunit-dev zlib1g-dev
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

If the compilation has run successfully, you should be able to run the `DEploidIBD` executable with:

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

*Recommended: Adding `DEploidIBD` to your bash PATH*
If we want to run the `dEploid` executable from anywhere, we need to add it to our `PATH` variable in bash. To do this, first make sure you are in the root directory of the `DEploidIBD` github repository. Then, run:
```
export PATH=$PATH:`pwd`
```
Now, you should be able to run `dEploid` from any folder.

*Additional resources: The `DEploid` manual page*
There are some additional instructions on how to run `DEploid` and the meaning of individual flags available in the github repository. First, navigate into the cloned DEploid directory on your local machine, e.g. `cd DEploid`. Then, run:
```
man docs/_build/man/dEploid.1
```

#### Step 4: Installing other tools required for the tutorial
In the upcoming sections of the tutorial we will make use of a command-line VCF manipulation tool called [bcftools](https://samtools.github.io/bcftools/bcftools.html). `bcftools` can be installed in a variety of ways, but I would recommend doing so using `conda`.

If you have not already, install `conda` following the appropriate instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Next, we will create a virtual environment for the tutorial that contains `bcftools` with the following command:
```
conda create -n deploid -c bioconda bcftools
```

For Mac you can also install `bcftools` by running:
```
brew install bcftools
```

Once this is done, you should be able to activate your environment with:
```
conda activate deploid
```
Running `bcftools` should display the help menu which includes a list of available subcommands. We are now ready to proceed with the first part of the tutorial.


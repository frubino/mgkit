.. _install-ref:

Installation
============

The library is both in `PyPI <https/www.pypi.org>`_ and `Bioconda <https://bioconda.github.io>`_. Any new version is first released on PyPI and then ported to the Bioconda repository. Installation is easier from Bioconda, since some dependencies (specifically `pysam <https://github.com/pysam-developers/pysam>`_) can have requirements that are difficult to resolve.

.. warning::

	The preferred version of Python to use is >=3.5 as this is the one I'm using MGKit with and Python 2.7 will not be mantained anymore starting with 1st January 2020. From MGKit 0.3.4 support for Python 3 was added and from 0.4.x you can expect Python 2.7 supoort to be gradually removed.

Using Conda
-----------

.. note::

	this was tested in conda version 4.5.4 and Python 3.6.5

Installing in a conda environment is relatively straight forward, and I recommend installing mgkit in its own environment::

	$ conda create --name mgkit mgkit

This will install `mgkit` in a virtual environment called `mgkit`. Conda will write on screen the way to activate the virtual environment, which may different according to your installation.

Docker Instance (with Jupyter Notebook)
---------------------------------------

A preconfigured Docker instance (user: mgkit, no password) has been configured at `Docker Hub <https://hub.docker.com/r/frubino/mgkit/>`_, including more packages for testing, available at Docker Hub (frubino/mgkit), with more instruction on its use available there. The version of MGKit targeted is the last development branch, but can be customised using the files available at `github <https://github.com/frubino/mgkit-docker-repo>`_, specifically in the `bootstrap.sh` file.

Pip and pipenv
--------------

`pipenv <https://pipenv.readthedocs.io/>`_ is the environment I use for developing and test, `conda` should be preferred for end users.

As of version 0.4.0, only python >=3.5 is used. The library has been ported to Python 3 (tested under version 3.5 and laters), but a layer of compatibility with Python 2.7 is used (the **future** package). Forward, version 0.3.4 is the last one that targets primarly Python 2.7 (but is partially compatible with Python 3). To test the version of Python installed use::

	$ python --version

.. warning::

	when installing pysam from source, the compilation of its libraries need system packages installed. Read the `INSTALL <https://github.com/pysam-developers/pysam/blob/master/htslib/INSTALL>`_ file to help you install the library

Assuming you have installed the required packages for `pysam` and other dependencies, `pipenv` can be used::

	$ pipenv install mgkit

and the environment activated with::

	$ pipenv shell

.. note::

	`enum34` is installed with Python version < 3.4

To recreate the `.c` files for the `mgkit.io.utils.sequences`, set the environment variable USE_CYTHON::

	$ export USE_CYTHON=TRUE

and then::

	$ python setup.py build

Using the repository
^^^^^^^^^^^^^^^^^^^^

The source code can also be obtained from the `Bitbucket repository <https://bitbucket.org/setsuna80/mgkit>`_.

Running Tests
---------------

Tests are performed with `tox <https://tox.readthedocs.io/en/latest/>`_, which you can use to run tests on specific versions of the interpreter. You can look at `tox.ini` in the distribution for the version tested.

The tests requires the `pytest` package and some plugins (`pytest-datadir` and `pytest-console-scripts`)::

	$ pip install pytest pytest-datadir pytest-console-scripts

You can run the tests with::

	$ python setup.py test

Some test won't be run if the required library/data is not found. Consult the output for more information.

Building Documentation
----------------------

The requirements are detailed in `docs_req.txt` of the repository. Other libraries:

* latex (for pdf output - `make latexpdf`)
* pandoc

Troubleshooting (for older versions)
------------------------------------

Some of the dependencies require available compilers to finish the installation. At the mimimum a system that provides the full GNU compiler suite, including a fortran compiler is required to install those dependencies by source.

If a compilation error is raised during installation, it's adviced to install each dependency manually. I'll try to keep this section updated, but there's not that many OS that I can keep working on (mostly MacOSX and GNU/Linux).

Older versions of MacOSX
^^^^^^^^^^^^^^^^^^^^^^^^

Version 10.19 of MacOSX comes with Python 2.7 installed. To install every dependency from source, however it's needed to install the *Xcode* app from the **App Store** which install the compilers, with the exception of `gfortran`. Another solution is using `Homebrew <http://brew.sh>`_ or `Macports <http://www.macports.org>`_, to install the compilers needed.

If you want to use Xcode, you need to install the gfortran compiler, with the package provided `here <http://gcc.gnu.org/wiki/GFortranBinariesMacOS>`_. This should be enough to install most packages from source.

.. warning::

	There seems to be a problem with `pandas` version 0.13.1 on MacOSX, with a segmentation fault happening when using DataFrames. The 0.14.1 version is the one tested.

.. note::

	if there's a problem building a python package because of a compile error, dealing with an unknown command line option, use::

		export ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future

	It's related to the clang toolchain included with Xcode

Matplotlib
**********

The tricky package to install in MacOSX is actually `matplotlib <http://matplotlib.org>`_, with one of many solutions being posted on `a disccusion on stackoverflow <http://stackoverflow.com/questions/4092994/unable-to-install-matplotlib-on-mac-os-x>`_. In our case, installing `freetype2` and `libpng` through Homebrew it's the less painful::

	$ brew install libpng freetype2

.. note::

	If you get a compilation error which refers to freetype2 in the `/opt/X11/` I found it easy to delete XQuartz installing matplotlib and then reinstall XQuartz.

	Or use::

		export PKG_CONFIG_PATH=/usr/local/Cellar/freetype/2.6_1/lib/pkgconfig/:/usr/local/Cellar/libpng/1.6.19/lib/pkgconfig/

	Note that the versions may be different.


Installing Scipy from source on Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A full description on how to install the scipy on Linux from source can be found at `this address <http://www.scipy.org/scipylib/building/linux.html>`_, be aware that the compilation of the `math-atlas` and `lapack` libraries takes a long time.

Installation in a virtual environment::

	# create virtual environment, if needed, otherwise activate the one desired
	virtualenv venv
	source venv/bin/activate
	# create temporary directory to compile math-atlas and lapack
	mkdir dep-build; cd dep-build
	wget http://www.netlib.org/lapack/lapack.tgz
	wget http://sourceforge.net/projects/math-atlas/files/Stable/3.10.2/atlas3.10.2.tar.bz2/download
	tar xfvj download
	cd ATLAS
	mkdir build; cd build
	../configure -Fa alg -fPIC --with-netlib-lapack-tarfile=../../lapack.tgz --prefix=$VIRTUAL_ENV
	make
	cd lib; make shared; make ptshared; cd ..
	make install

This will compile math-atlas with full lapack support in the virtual environment; change the `--prefix=$VIRTUAL_ENV` to `--prefix=$HOME` if you want to install the dependencies in you home directory.

Notes
-----

Not all packages are required to use the part of the library, but it's
recommended to install all of them. Requirements are bound to change, but pandas, scipy,
numpy, pysam and matplotlib are the bases of the library.

To avoid problems with the system installation, I suggest using the excellent
`virtualenv <http://www.virtualenv.org/>`_. This will avoid problems with
installing packages system-wide and breaking a working installation.


.. rubric:: Footnotes

.. [#] http://www.pip-installer.org/en/latest/user_guide.html#user-installs

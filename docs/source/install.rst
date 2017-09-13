.. _install-ref:

Installation
============

Docker Instance (with Jupyter Notebook)
---------------------------------------

A preconfigured Docker instance (user: mgkit, no password) has been configured using the instructions in :ref:`install-ubuntu`, including more packages for testing, available at Docker Hub (frubino/mgkit)::

	$ docker run -p 8881:8888 -v host-dir:/home/mgkit/notebooks/ -it frubino/mgkit

This command (assuming that Docker is already installed), will pull the instance and present with a bash terminal. IPython and Jupyter are installed and can be used. For a recap of the options:

* `-p 8888:8888` instruct to open the port 8888 on the host
* `-v` mount the directory `host-dir` into the virtual machine */home/mgkit/notebooks/* directory
* `-it` opens an interactive shell

The port to open is port 8881 on the host for using the Jupyter Notebook. The port can be change to what best fits the user. The same applies to the working directory.

After entering the virtual machine prompt, the command that is to be used to launch the notebook is::

	$ jupyter notebook --ip=0.0.0.0

After which in a browser it will possible to access it at *localhost:8881*.

.. note::

	The Dockerfile used to build the instance is available in the repository at: docs/source/extra/Dockerfile


Requirements
------------

The library is written completely in `Python <http://www.python.org>`_ and has been tested with both version 2.6.x and 2.7.x. It has not been tested with Python 3.x, but the authors are trying to write code that is potentially runnable with the use of the `2to3` tool.

Most UNIX systems provide a version of Python installed. Latest versions of MacOSX provides Python 2.7.x, while most Linux variants may have different versions installed, sometimes version 3.x, but they usually provide python 2.7.x packages (`Archlinux <https://www.archlinux.org/>`_ uses version 3.x by default, but a `python2` package is provided). You can check the Python version installed with::

	$ python --version

The library requires these Python packages:

The installation dependencies are flexible, with only *numpy* (version 1.9.2 or later) as being **required**::

	$ pip install mgkit

To install every needed packages, you can use::

	$ pip install mgkit[full]

Other options can be found in the *setup.py* file, in the *extra_requires* dictionary.

The optional dependencies includes:

* scipy
* pandas
* matplotlib
* pysam
* argparse (if Python 2.6 is installed, part of the standard library from 2.7)
* joblib: for script `translate_seq`, to use multiple processors
* HTSeq
* semidbm
* pymongo
* `goatools <https://github.com/tanghaibao/goatools>`_ (required by :mod:`mgkit.mappings.go`), has package `fisher` as a dependency
* rpy2 >= 2.3.8 (required by :mod:`mgkit.utils.r_func`)

.. _install-ubuntu:

Installing on Ubuntu Server 16.04
---------------------------------

You'll need to install the following packages with `apt-get`::

	$ apt-get install -y velvet bowtie2 python-pip python \
	  virtualenv python-dev zlib1g-dev libblas-dev \
	  liblapack-dev gfortran libfreetype6-dev libpng-dev \
	  fontconfig pkg-config

Create a virtual environment to ensure that the correct library versions are installed as explained in :ref:`install-virtualenv`.

Using pip
---------

All dependencies are usually installed either through a package system provided by the running OS or through the `pip <http://www.pip-installer.org/>`_ installer. If you're using a system that's shared with other people, you may not be able to install the dependencies system-wide, in which case the `--user` option in `pip` may solve the problem [#]_.

A system-wide installation with `pip` can be done with::

	$ pip install path/to/library

while a user install is done with::

	$ pip install --user path/to/library

all requirements will be downloaded/installed.

.. _install-virtualenv:

Using virtualenv
^^^^^^^^^^^^^^^^

`virtualenv <http://www.virtualenv.org/>`_ is a system that is used to isolate a Python installation, to make sure no conflicts arise with multiple packages. It's handy if you're developing or testing an application/library, as it provides a clean environment.

Assuming you've already installed `virtualenv`, a virtual environment can be created with::

	$ virtualenv -p python2 mgkit-env

which creates a virtual environment in `mgkit-env`, with the interpreter used being the one linked to `python2`. Activating the environment requires using::

	$ source mgkit-env/bin/activate

assuming you're in the same directory where you created the environment. The pip packager is installed by default with it, so we're going to use it to install the library if you have downloaded it already::

	$ pip install path/to/library

or getting the last version from `PyPI <https://pypi.python.org/pypi>`_::

	$ pip install mgkit

You can also install a specific version::

	$ pip install mgkit==0.2.0

Using the repository
^^^^^^^^^^^^^^^^^^^^

The source code can also be obtained from the `Bitbucket repository <https://bitbucket.org/setsuna80/mgkit>`_.

Running Tests
---------------

The tests requires the `nosetests` package::

	$ pip install nose

and the package `yanc` is used for coloring the output. If you don't want to install it you can edit the `setup.cfg` and `setup.py` files in the source distribution and delete the `with-yanc` before running the tests.

You can run the tests with::

	$ python setup.py nosetests

Some test won't be run if the required library/data is not found. Consult the output for more information.

Building Documentation
----------------------

Needs sphinx >=1.2.2

* sphinx_rtd_theme
* actdiag
* sphinxcontrib-actdiag
* blockdiag
* sphinxcontrib-blockdiag
* sphinxcontrib-napoleon (we'll be part of sphinx 1.3, needed until then)
* sphinx-argparse

Other libraries:

* graphviz
* latex (for pdf output - `make latexpdf`)

Troubleshooting
---------------

Some of the dependencies require available compilers to finish the installation. At the mimimum a system that provides the full GNU compiler suite, including a fortran compiler is required to install those dependencies by source.

If a compilation error is raised during installation, it's adviced to install each dependency manually.

I'll try to keep this section updated, but there's not that many OS that I can keep working on (mostly MacOSX and GNU/Linux).

HTSeq
^^^^^

Sometimes HTSeq or numpy fails to install in a clean environment; it's advised to install numpy first::

	$ pip install numpy 

and then reissue the library installation::

	$ pip install path/to/library

MacOSX
^^^^^^

The version of MacOSX is 10.9 that comes with Python 2.7 installed. To install every dependency from source, however it's needed to install the *Xcode* app from the **App Store** which install the compilers, with the exception of `gfortran`. Another solution is using `Homebrew <http://brew.sh>`_ or `Macports <http://www.macports.org>`_, to install the compilers needed.

If you want to use Xcode, you need to install the gfortran compiler, with the package provided `here <http://gcc.gnu.org/wiki/GFortranBinariesMacOS>`_. This should be enough to install most packages from source.

.. warning::

	There seems to be a problem with `pandas` version 0.13.1 on MacOSX, with a segmentation fault happening when using DataFrames. The 0.14.1 version is the one tested.

.. note::

	if there's a problem building a python package because of a compile error, dealing with an unknown command line option, use::

		export ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future

	It's related to the clang toolchain included with Xcode

Scipy
*****

There are different solutions available if you have trouble installing the dependencies on MacOSX, one of which is hosted `on this page <http://fonnesbeck.github.io/ScipySuperpack/>`_, but installing from source is another option, provided that the Xcode and gfortran are installed.

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

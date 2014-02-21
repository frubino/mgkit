.. _install-ref:

Installation
============

Requirements
------------

The library is written completely in `Python <http://www.python.org>`_ and has been tested with both version 2.6.x and 2.7.x. It has not been tested with Python 3.x, but the authors are trying to write code that is potentially runnable with the use of the `2to3` tool.

Most UNIX systems provide a version of Python installed. Latest versions of MacOSX provides Python 2.7.x, while most Linux variants may have different versions installed, sometimes version 3.x, but they usually provide python 2.7.x packages (`Archlinux <https://www.archlinux.org/>`_ uses version 3.x by default, by a `python2` package is provided). You can check the Python version installed with::
	
	$ python --version

The library requires these Python packages:

* numpy >= 1.7.1
* scipy >= 0.13.0
* pandas >= 0.12.0
* HTSeq >= 0.5.4p5
* matplotlib >= 1.3.1
* pysam >= 0.7.7 (required by :mod:`mgkit.align`)
* goatools (required by :mod:`mgkit.mappings.go`)
* rpy2 >= 2.3.8 (required by :mod:`mgkit.utils.r_func`)
* argparse (if Python 2.6 is installed, part of the standard library from 2.7)
 
Using pip
---------

All dependencies ca usually installed either through a package system provided by the running OS or through the `pip <http://www.pip-installer.org/>`_ installer. If you're using a system that's shared with other people, you may not be able to install the dependencies system-wide, in which case the `--user` option in `pip` may solve the problem [#]_.

A system-wide installation with `pip` can be done with::

	$ pip install path/to/library

while a user install is donw with::

	$ pip install --user path/to/library

all requirements we'll be downloaded/installed.

Using virtualenv
----------------

`virtualenv <http://www.virtualenv.org/>`_ is a system that is used to isolate a Python installation, to make sure no conflicts arise with multiple packages. It's handy if you're developing or testing an application/library, as it provides a clean environment. 

Assuming you've already installed `virtualenv`, a virtual environment can be created with::

	$ virtualenv -p python2 mgkit-env

which create a virtual environment in `mgkit-env`, with the interpreter used being the one linked to `python2`. Activating the environment requires using::

	$ source mgkit-env/bin/activate

assuming you're in the same directory where you created the environment. The pip packager is installed by default with it, so we're going to use it to install the library::
	
	$ pip install path/to/library

Using the repository
^^^^^^^^^^^^^^^^^^^^

bitbucket link

Runnining Tests
---------------

The tests requires the `nosetests` package::

	$ pip install nosetests

and the package `yanc` is used for coloring the output. If you don't want to install it you can edit the `setup.cfg` and `setup.py` files in the source distribution and delete the `with-yanc` before running the tests.

You can run the tests with::

	$ python setup.py nosetests

Some test won't be run if the required library/data is not found. Consult the output for more information.

Building Documentation
----------------------

Needs sphinx >=1.2, sphinxcontrib.napoleon (if sphinx <1.3), rst2pdf, sphinxcontrib.blockdiag, sphinxcontrib.actdiag, sphinx.ext.graphviz, plus latex installed (needs changing the conf.py file)

Troubleshooting
---------------

Some of the dependencies requires available compilers to finish the installation. At the mimimum a system that provides the full gnu compiler suite, including a fortran compiler is required to install those dependencies by source.

If a compilation errors is raised during installation, it's adviced to install each dependency by hand.

I'll try to keep this section updated, but there's not that many OS that I can keep working on (mostly MacOSX and Linux).

HTSeq
^^^^^

Sometimes HTSeq or numpy fails to install in a clean environment; it's advised to install numpy first::

	$ pip install numpy 

and then reissue the library installation::

	$ pip install path/to/library

MacOSX
^^^^^^

There are different solutions available if you have trouble installing the dependencies on MacOSX, one of which is hosted `on this page <http://fonnesbeck.github.io/ScipySuperpack/>`_,

The version tested is MacOSX (10.9) that comes with Python 2.7 installed. To install every dependency from source, however it's needed to install the *Xcode* app from the app store which install the compilers, with the exception of `gfortran` or a solution line `Homebrew <http://brew.sh>`_ or `Macports <http://www.macports.org>`_, which can be used to install the compilers needed.

Either solution is fine, but if you wnat to use Xcode, you need to install the gfortran compiler, with the package provided `here <http://gcc.gnu.org/wiki/GFortranBinariesMacOS>`_. This should be fine to install most packages from source.

The tricky package to install in MacOSX is actually `matplotlib <http://matplotlib.org>`_, with one of many solutions being posted on `a disccusion on stackoverflow <http://stackoverflow.com/questions/4092994/unable-to-install-matplotlib-on-mac-os-x>`_. In our case, installing `freetype2` and `libpng` through Homebrew it's the less painful::

	$ brew install libpng freetype2

.. note::

	If you get a compilation error which refers to freetype2 in the `/opt/X11/` I found it easy to delete XQuartz installing matplotlib and then reinstall XQuartz.

Installing scipy on Linux
^^^^^^^^^^^^^^^^^^^^^^^^^

In case you can't install scipy from the source, because of a compile error, you may try the solution on `stackoverflow  <http://stackoverflow.com/questions/7496547/python-scipy-needs-blas>`_

Remember to add **-fPIC** to the compilation options in LAPACK to the **make.inc**.

Notes
-----

Not all packages are required to use the part of the library but it's
recommended to do it. Requirements are bound to change, but pandas, scipy,
numpy, pysam and matplotlib are the bases of the library.

To avoid problems with the system installation, I suggest using the excellent
`virtualenv <http://www.virtualenv.org/>`_. This will avoid problems with
installing packages system-wide and breaking a working installation.


.. rubric:: Footnotea

.. [#] http://www.pip-installer.org/en/latest/user_guide.html#user-installs

.. _sec-installation:

============
Installation
============

``dEploid`` is written in C++. We use `autotools` for compiling and installing.

**************
Stable Release
**************

The latest version of ``dEploid`` is available at `github <https://github.com/mcveanlab/DEploid/releases/latest>`_.


*******************************
Development Version From GitHub
*******************************

You can also install ``dEploid`` directly from the git repository (`tar <https://api.github.com/repos/mcveanlab/DEploid/tarball>`_, `zip <https://api.github.com/repos/mcveanlab/DEploid/zipball>`_). Here, you will need ``autoconf``, check whether this is already installed by running:

.. code-block:: bash

    $ which autoconf

On Debian/Ubuntu based systems:

.. code-block:: bash

    $ apt-get install build-essential autoconf autoconf-archive libcppunit-dev


On Mac OS:

.. code-block:: bash

    $ port install automake autoconf autoconf-archive cppunit


Afterwards you can clone the code from the github repository,

.. code-block:: bash

    $ git clone git@github.com:mcveanlab/DEploid.git
    $ cd DEploid

and build the binary using

.. code-block:: bash

    $ ./bootstrap
    $ make

or install with ``make install``.

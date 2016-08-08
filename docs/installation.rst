.. _sec-installation:

============
Installation
============

``DEploid`` is written in C++.

**************
Stable Release
**************

.. todo::
    You can download the latest stable release packaged for a variety of different platform from


*******************************
Development Version From GitHub
*******************************

You can also install ``DEploid`` directly from the git repository. Here, you need to install ``autoconf`` first:

On Debian/Ubuntu based systems:

.. code-block:: bash

    apt-get install build-essential autoconf autoconf-archive libcppunit-dev


On Mac OS:

.. code-block:: bash

    port install automake autoconf autoconf-archive cppunit


Afterwards you can build the binary using

.. code-block:: bash

    ./bootstrap
    make

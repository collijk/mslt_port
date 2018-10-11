.. _installation:

Installation
============

Python
------

Before you can install Vivarium_ or the MSLT_ framework, you will need a
working installation of Python.
If you don't have Python installed the best option is to use Anaconda_, which
is available for Windows, Mac OS X, and Linux, and includes all of the
necessary scientific packages.

Vivarium
--------

Once you have a working Python installation, you can install Vivarium_:

.. code-block:: shell

   pip install vivarium

You can verify that Vivarium_ was successfully installed by opening the Python
interpreter and running:

.. code-block:: python

   import vivarium
   vivarium.__version__

This should display the version of Vivarium_ that you have installed.

MSLT framework
--------------

You can then install the MSLT_ framework by making a local copy of the
repository, opening a command prompt, and running:

.. code-block:: shell

   cd path/to/mslt_port
   pip install .

You can verify that the MSLT_ framework was successfully installed by opening
the Python interpreter and running:

.. code-block:: python

   import mslt_port

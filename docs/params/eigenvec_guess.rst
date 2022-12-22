*eigenvec_guess*
======================

Syntax
""""""

.. parsed-literal::

   eigenvec_guess = arg

* arg = character


Default
"""""""

.. code-block:: fortran

   eigenvec_guess = " "


Description
"""""""""""

Filename where the eigenvector guess is read. This option can be used when a specific eigenvector should be used to start the calculation, for example to refine a saddle point when the eignevector is known.

The file format is *xyz*, vector read from file is used as-is.


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

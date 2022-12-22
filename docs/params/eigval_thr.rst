*eigval_thr*
======================

Syntax
""""""

.. parsed-literal::

   eigval_thr = arg

* arg = real


Default
"""""""

.. code-block:: bash

   eigval_thr = -0.01


Description
"""""""""""

Threshold for eigenvalue, which determines when to start following the eigenvector.
It specifies the line of *inflection* of the curvature of the PES, the region below it can be considered as the basin around some minimum.

It is generally a good idea to specify a small, but not-too-small negative value, in order to be sure to exit the basin region of the PES before starting to follow the eigenvector direction.

In units of ``engine_units``.


Unexpected behavior
"""""""""""""""""""

Can vary greatly from system to system, and precision of the calculated eigenvalue.


Related commands
""""""""""""""""

:doc:`nsmooth`

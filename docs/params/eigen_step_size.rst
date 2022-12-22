*eigen_step_size*
======================

Syntax
""""""

.. parsed-literal::

   eigen_step_size = arg

* arg = real


Default
"""""""

.. code-block:: bash

   eigen_step_size = 0.2


Description
"""""""""""

The limit to the maximum size of the displacement with eigenvector.

Actual step size is calculated from the force and eigenvalue, see :file:`/src/artn.f90`, variable ``current_step_size``.


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

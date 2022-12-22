*push_step_size*
======================

Syntax
""""""

.. parsed-literal::

   push_step_size = arg

* arg = real


Default
"""""""

.. code-block:: fortran

   push_step_size = 0.4


Description
"""""""""""

Maximum size of a component in the initial displacement vector.

This step size is kept fixed as long as the push direction follows the initial push vector. When pushing with eigenvector, the step size is regulated.


Unexpected behavior
"""""""""""""""""""

Can be limited by the E/F engine parameters. For example, when using a large ``push_step_size`` with LAMMPS, the parameter ``dmax`` of FIRE in LAMMPS should be accordingly modified when invoking ``fix artn`` in the LAMMPS input, *e.g.*:

.. code-block:: bash

   fix 3 all artn dmax 8.0



Related commands
""""""""""""""""

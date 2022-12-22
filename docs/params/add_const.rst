*add_const*
======================

Syntax
""""""

.. parsed-literal::

   add_const = arg

* arg = real(4,:) array


Default
"""""""

.. code-block:: fortran

   add_const(:,:) = 0.0


Examples
""""""""

.. code-block:: bash

   1.
   push_ids = 5, 2
   add_const(:, 2) = 1.0, 0.0, -1.0, 30.0

   2.
   push_ids = 6, 13
   add_const(:, 13) = 0.5, -0.5, 0.5, 0.0
   add_const(:, 6) = 0.3, 1.0, -1.0, 90.0



Description
"""""""""""

Constraint on the initial push on each atom: 3-component vector, and 1 solid angle in degrees.
The 3-component vector gets normalized, and does not have any effect on the ``push_size``.

Explanation of the examples above:
 1. generate push on atom index 2, which is a random vector within the cone specified by the vector (1, 0, -1), and the angle 30 degrees from its apex, where the apex is the atomic position of atom index 2. The initial push is also generated for atom 5, but without contraints.
 2. push on atoms 13 and 6, atom 13 has 0.0 angle, the push vector is exact (0.5, -0.5, 0.5); atom 6 has angle 90 degrees, which is any vector from the hemisphere around (0.3, 1.0, -1.0).


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

:doc:`push_ids`, :doc:`push_mode`

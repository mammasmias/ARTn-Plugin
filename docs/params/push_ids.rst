*push_ids*
======================

Syntax
""""""

.. parsed-literal::

   push_ids = arg

* arg = integer(:) array


Default
"""""""

.. code-block:: fortran

   push_ids(:) = 0


Description
"""""""""""

List of atom indices with nonzero components in the initial push vector.
This array is used when ``push_mode = 'list'``, and ``push_mode = 'rad'``.

In order to constrain the randomly generated initial push vector in some direction, use the parameter ``add_const``.


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

:doc:`push_mode`, :doc:`add_const`

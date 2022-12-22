*dist_thr*
======================

Syntax
""""""

.. parsed-literal::

   dist_thr = arg

* arg = real


Default
"""""""

.. code-block:: bash

   dist_thr = 0.0


Example
"""""""

.. code-block:: bash

   push_mode = 'rad'
   push_ids = 7
   add_const(:,7) = 0.0, 1.0, 1.0, 15.0
   dist_thr = 3.2



Description
"""""""""""

Generate push on all atoms within the radius from an atom in ``push_ids``, used in combination with ``push_mode = 'rad'``.


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

:doc:`push_mode`, :doc:`push_ids`, :doc:`add_const`

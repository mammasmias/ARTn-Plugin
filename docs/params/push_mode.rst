*push_mode*
======================

Syntax
""""""

.. parsed-literal::

   push_mode = arg

* arg = character


Default
"""""""

.. code-block:: bash

   push_mode = 'all'


Description
"""""""""""

Specify atomic indices on which to generate a non-zero initial push vector.
Possible values are ``'all'``, ``'list'``, ``'rad'``, or ``'file'``.

When using ``push_mode='list'``, the initial push is generated only on the atoms present in the list of atomic indices, which needs to be specified in ``push_ids``.

When using ``push_mode='rad'``, the initial push is generated on all atoms within the radial distance of atoms in the list of indices ``push_ids``, the radial cutoff needs to be provided in ``dist_thr``.

When using ``push_mode='file'``, the initial push is read from a file, whose filename needs to be provided in ``push_guess``.

Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

:doc:`push_ids`, :doc:`dist_thr`, :doc:`push_guess`.

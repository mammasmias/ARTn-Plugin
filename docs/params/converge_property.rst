*converge_property*
======================

Syntax
""""""

.. parsed-literal::

   converge_property = arg

* arg = character, possible values ``maxval``, ``norm``.


Default
"""""""

.. code-block:: bash

   converge_property = "maxval"


Description
"""""""""""

Specify how to test convergence of the forces:
 - ``converge_property = 'maxval'`` the convergence will be tested by ``MAXVAL( ABS( force ) )``;
 - ``converge_property = 'norm'`` the convergence will be tested by ``NORM2( force )``.


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

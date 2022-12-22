*engine_units*
======================

Syntax
""""""

.. parsed-literal::

   engine_units = arg/mode

* arg = character, possible values ``'qe'``, ``'lammps'``
* mode = character, when arg=``'lammps'`` possible values are ``'metal'``, ``'real'``, ``'lj'``


Default
"""""""

.. code-block:: bash

   engine_units = 'qe'


Examples
""""""""

.. code-block:: bash

   engine_units = 'qe'
   engine_units = 'lammps/metal'
   engine_units = 'lammps/real'


Description
"""""""""""

Select the units according to the E/F engine of choice:
 - ``'qe'`` units used in Quantum ESPRESSO: Rydberg, bohr, a.u.time
 - ``'lammps/mode'`` units defined in LAMMPS depending on ``'mode'``.

In order to use E/F engine units which are not listed here, you need to implement them on your own. See file :file:`/src/units.f90`.


Unexpected behavior
"""""""""""""""""""


Related commands
""""""""""""""""

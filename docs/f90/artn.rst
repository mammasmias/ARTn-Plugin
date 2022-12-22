.. _f90_artn:

artn
=======

This is manually added Fortran-style block, that gets interpreted by ``sphinx-fortran`` module.
For more info see `this link <https://sphinx-fortran.readthedocs.io/en/latest/user.domain.html#>`_.

.. f:subroutine:: artn(force, nat, fout)

   description of subroutine in .rst file written manually, not captured from source code.

   :p float force(3*nat) [in]: engine force
   :p integer nat [in]: number of atoms
   :r float fout(3*nat): modified force

   :calledfrom: :file:`artn_QE.f90`
   :callto: :file:`lanczos.f90`

--------------

The following is captured from Doxygen and merged with the ``breathe`` software.
For documentation see `this website <https://breathe.readthedocs.io/en/latest/index.html>`_.

.. doxygenfile:: artn.f90
   :project: plugin-ARTn

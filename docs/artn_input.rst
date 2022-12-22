List of ARTn input parameters
=============================

The ARTn parameters are given in the file `artn.in`. That file is formatted as FORTRAN NAMELIST, called `ARTN_PARAMETERS` as for example:

.. code-block:: bash

   &ARTN_PARAMETERS
                ... specify parameters ...
   /


All parameters available in pARTn are listed below, grouped by the part of ARTn algorithm they affect.

.. toctree::
   :maxdepth: 1
   :caption: I/O control and options

   params/verbose
   params/engine_units
   params/struc_format_out
   params/lrestart
   params/lpush_final
   params/lmove_nextmin
   params/zseed


.. toctree::
   :maxdepth: 1
   :caption: Control initial push

   params/push_mode
   params/push_ids
   params/add_const
   params/dist_thr
   params/push_step_size
   params/push_guess
   params/ninit


.. toctree::
   :maxdepth: 1
   :caption: Control the Lanczos algorithm

   params/lanczos_max_size
   params/lanczos_disp
   params/lanczos_eval_conv_thr


.. toctree::
   :maxdepth: 1
   :caption: Control the eigenvector push

   params/eigval_thr
   params/eigen_step_size
   params/eigenvec_guess
   params/nsmooth
   params/neigen
   params/nnewchance


.. toctree::
   :maxdepth: 1
   :caption: Control the perpendicular relaxation

   params/nperp
   params/nperp_limitation


.. toctree::
   :maxdepth: 1
   :caption: Control convergence

   params/forc_thr
   params/converge_property





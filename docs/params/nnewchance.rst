*nnewchance*
======================

Syntax
""""""

.. parsed-literal::

   nnewchance = arg

* arg = integer


Default
"""""""

.. code-block:: bash

   nnewchance = 0


Description
"""""""""""

Number of times a research is allowed to re-initialize the first Lanczos vector.

When the system has reached an area on the PES where there are positive curvatures but no minimum, as an :math:`x^3` curve. It is not clear for the moment when this happens and how it can be detected or avoided.

Unexpected behavior
"""""""""""""""""""

Can lead to unconnected paths.


Related commands
""""""""""""""""

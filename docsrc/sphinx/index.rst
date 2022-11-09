
RePrimAnd Documentation
=======================

Recovery of Primitives And EOS framework


Overview
--------

`RePrimAnd` is a support library for numerical simulations of general 
relativistic magnetohydrodynamics. 
It provides methods for recovering primitive variables like pressure
and velocity from the variables evolved in quasi-conservative 
formulations. Further, it provides a general framework for handling
matter equations of state, a TOV solver, and tools for precomputed
NS sequences.

A description of the recovery algorithm and tests of the library is given in 
the accompanying article :footcite:p:`Kastaun2021`.

The latest public versions can be found on the 
`github page <https://github.com/wokast/RePrimAnd>`_.
Releases are also archived on Zenodo :footcite:p:`Zenodo:Reprimand`.

Developers
^^^^^^^^^^
`RePrimAnd` was written by Dr. Wolfgang Kastaun 
<physik@fangwolg.de>.
Thanks for additional testing go to Jay Vijay Kalinani.


License
^^^^^^^

.. raw:: html

   <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

Please cite the PRD article :footcite:p:`Kastaun2021`.
and the DOI on Zenodo :footcite:p:`Zenodo:Reprimand`
when publishing results obtained using the library.


References
^^^^^^^^^^

.. footbibliography::



Contents
--------

.. toctree::
   :maxdepth: 2

   installing
   notation
   c2p_imhd
   eos_thermal
   eos_barotr
   tov_solver
   ns_seqs
   little_helpers
   python_interf

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

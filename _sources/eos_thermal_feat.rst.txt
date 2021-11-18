Features
========

The library contains different types of 3-parametric EOS that can 
be accessed through a generic :ref:`interface<thermal_interface>`. 
It is quite simple to add custom EOS without modifying the library.

The interface provides the following functions:

* :math:`P(\rho,\epsilon,Y_e)`
* :math:`c_s(\rho,\epsilon,Y_e)`
* :math:`\frac{\partial P}{\partial \rho}(\rho,\epsilon,Y_e)`
* :math:`\frac{\partial P}{\partial \epsilon}(\rho,\epsilon,Y_e)`
* :math:`s(\rho,\epsilon,Y_e)` (optional)
* :math:`T(\rho,\epsilon,Y_e)` (optional)

Optionally, the above is also available in terms of temperature 
:math:`T` instead of :math:`\epsilon`. 

Further, each EOS comes with a validity region of the form

.. math::

   \rho_\mathrm{min} &\le \rho \le \rho_\mathrm{max}   \\
   Y_{e,\mathrm{min}} &\le Y_e \le Y_{e,\mathrm{max}}      \\
   \epsilon_\mathrm{min}(\rho, Y_e) &\le \epsilon 
     \le \epsilon_\mathrm{max}(\rho, Y_e)  \\
   T_\mathrm{min}(\rho, Y_e) &\le T  \le T_\mathrm{max}(\rho, Y_e) \\

In addition, the global minimum :math:`h_0` of the enthalpy over the 
valid region is given for each EOS. 

   
The framework guarantees some physical constraints:

.. math::
   
   0 &\le c_s < 1 \\
   0 &\le P \le \rho_E  \\
   \rho_E &\ge 0 \\
   h_0 &> 0
   


.. note::

   Only :math:`P(\rho,\epsilon,Y_e), h_0`, and the validity region
   are used for the primitive recovery, the other functions are
   intended for use in numerical relativity evolution codes.


    



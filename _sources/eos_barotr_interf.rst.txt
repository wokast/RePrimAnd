.. _barotr_interface:

Interface
---------


The barotropic EOS are represented
by objects of type :cpp:class:`~EOS_Toolkit::eos_barotr`. 
These objects 
represent the generic EOS interface and do not implement any EOS. 
Instead, they forward any request to an actual implementation of a 
particular EOS. To obtain an 
:cpp:class:`~EOS_Toolkit::eos_barotr`
object representing a particular EOS, one uses functions such as 
:cpp:func:`~EOS_Toolkit::make_eos_barotr_poly` 
provided for each type of EOS.


Calling the EOS is a two stage process: first, one needs to specify
the state of matter. This can be done either in terms of mass density,
or in terms of the pseudo enthalpy :math:`g`.
For this, use the methods 
:cpp:func:`~eos_barotr::at_rho` and
:cpp:func:`~eos_barotr::at_gm1`. 
Both return an object of type
:cpp:class:`~EOS_Toolkit::eos_barotr::state`. 
All EOS always provide both forms.

Second, one can use methods of this object to obtain mass density, 
pseudo enthalpy, pressure, 
soundspeed, internal specific energy, temperature, and specific entropy. 
Further, one can use the state object like a boolean 
variable to check if the state is in the valid regime of the EOS. 
Trying to compute any of the above from an invalid state will throw 
an exception. The following snipped demonstrates the EOS use:

.. code-block:: cpp

   auto s = eos.at_rho(rho);
   if (s) {
     auto pressure   = s.press()
     auto soundspeed = s.csnd()
   }
   
   
The functionality above is also provided by another set of EOS methods 
with names such as :cpp:func:`~eos_barotr::press_at_rho`, 
:cpp:func:`~eos_barotr::press_at_gm1`, etc. Those do not throw 
exceptions when called outside the valid ranges, but instead return
NANs. When computing several quantities for the same state, the 
alternative syntax might be less efficient because the validity has to 
be checked each time instead of only once as for the primary syntax.

The :cpp:class:`~EOS_Toolkit::eos_barotr` class
also provides methods to query the valid ranges of mass density and 
pseudo enthalpy. For convenience, there are methods to check if 
given parameter is inside the valid range.

.. tip::

   The interface objects are designed to be used as ordinary variables,
   which can be copied around cheaply and without worrying about ownership 
   of any resource allocations. Internally, they just contain a shared 
   pointer to an actual (immutable) implementation.

.. note::

   The EOS interface is thread-safe in the sense that the same EOS 
   oject can be used in multiple threads. Since the EOS implementations
   are immutable, data-races cannot occur. Copying an EOS object inside
   a thread is save but might affect performance since the copying is 
   an atomic operation. Passing the EOS inside an innermost loop should
   better be done by reference instead of value. 
   
   


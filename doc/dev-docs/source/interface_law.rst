Interface Law
=============

cRacklet has several interfacial behavior implemented which are described in this section. These law are listed in :cpp:enum:`FractureLawType <FractureLawType>`.

Mixed-mode cohesive law
-----------------------

This cohesive law has its strength decreasing linearly with the opening. The opening considered here is the norm of the opening displacement, thus taking into account both normal and shear componenets. The critical displacement :math:`\delta_c` is the critical displacement required to transition from the peak stress :math:`\tau_c` to the residual value :math:`\tau_r` (By default is 0). The shear and normal components can be prescribed independantly. While the norm of the opening :math:`||delta||` is lower than :math:`\delta_c`, the strength is given by:  

.. math::
   \tau = \tau_c \left(1- ||\delta||/\delta_c \right)


Rate and state
--------------

The rate and state framework involve two functionnals: one for the friction coefficient :math:`f(v,\phi)` and one evolution law for the state variable :math:`\dot\phi = g(v,\phi)`. Several variations of the friction law and the evolution law are available in cRacklet. Their exact formulations are summarized here: 

Pure Weakening Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^

This formulation is based on the one originaly proposed by `Dieterich (1979) <https://doi.org/10.1029/JB084iB05p02161>`_ and `Ruina (1983) <https://doi.org/10.5194/npg-15-1-2008>`_ :

.. math::
   f(v,\phi) = f_0 + a \log \left(v/v_* \right) + b f_0 \log \left(\phi / \phi_* \right)

Standard Formulation
^^^^^^^^^^^^^^^^^^^^

This formulation is a modified version of the original formulation proposed by `Dieterich (1979) <https://doi.org/10.1029/JB084iB05p02161>`_ and `Ruina (1983) <https://doi.org/10.5194/npg-15-1-2008>`_ , with the addition of :math:`+1` in the logarithm.

.. math::
   f(v,\phi) = f_0 + a \log \left(1+v/v_* \right) + b f_0 \log \left(1 + \phi / \phi_* \right)

Regularized Formulation
^^^^^^^^^^^^^^^^^^^^^^^

This formulation is a generic N-shape friction law, as supported by experimental observations `Bar-Sinai & al. (2014) <https://doi.org/10.1002/2013JB010586>`_ . See `Brener & al. (2018) <https://doi.org/10.1103/PhysRevLett.121.234302>`_ for a discussion of the physical sense of each term and the comparison with more conventional rate and state friction law.. The steady-state friction is strenghtening at low and high velocities, and has a velocity-weakening branch at intermediate velocities.

.. math::
   f(v,\phi) = \left( 1 + b \log \left(1+ \frac{\phi}{\phi_*} \right) \right) \left( \frac{f_0}{  \sqrt{\left(  1+v_0^2 / v^2 \right)} } + a \log \left( 1+\frac{v}{v_*} \right) \right)

Regularized Weakening
^^^^^^^^^^^^^^^^^^^^^

This formulation is derived from the N-shape one, with the omission of the :math:`+1` term in the :math:`log` of :math:`\phi` term, resulting in the absence of the velocity-strengthening branch at high velocites.
   
.. math::
   f(v,\phi) = \left( 1 + b \log \left(\frac{\phi}{\phi_*} \right) \right) \left( \frac{f_0}{  \sqrt{\left(  1+v_0^2 / v^2 \right)} } + a \log \left( 1+\frac{v}{v_*} \right) \right)

Aging Law
^^^^^^^^^

The original aging law proposed by `Dieterich (1979) <https://doi.org/10.1029/JB084iB05p02161>`_.

.. math::
   g(\phi) = 1 - \frac{v \phi}{D}

Regularized Aging Law
^^^^^^^^^^^^^^^^^^^^^

A regularization of the slip law, ensuring that for vanishingly small steady-state velocity the value of the state :math:`\phi` saturates to a finite value :math:`D / v_*` after long times instead of diverging.

.. math::
   g(\phi) = 1 - \frac{v \phi}{D} \sqrt{1 + \left(\frac{v_*}{v}\right)^2}

Slip Law
^^^^^^^^

The slip law proposed by `Ruina (1983) <https://doi.org/10.5194/npg-15-1-2008>`_.

.. math::
   g(\phi) = - \frac{v \phi}{D} \log \left( \frac{v \phi}{D} \right)

:cpp:enum:`FractureLawType <_rate_and_state>`
   
Friction law
------------

Interface Law
=============

cRacklet has several interfacial behavior implemented which are described in this section. These law are listed in :cpp:enum:`FractureLawType <FractureLawType>`.

Mixed-mode cohesive law
-----------------------

This cohesive law has its strength decreasing linearly with the opening. The opening considered here is the norm of the opening displacement, thus taking into account both normal and shear componenets. The critical displacement :math:`\delta_c` is the critical displacement required to transition from the peak stress :math:`\tau_c` to the residual value :math:`\tau_r` (By default is 0). The shear and normal components can be prescribed independantly. While the norm of the opening :math:`||\delta||` is lower than :math:`\delta_c`, the strength is given by:  

.. math::
   \tau^{str} = \tau_c \left(1- ||\delta||/\delta_c \right)

Friction law for cohesive law
-----------------------------

In addition to the mixed-mode cohesive law, one can defined a friction law to handle cases when the surfaces are in contact with each others. 

Coulomb friction
^^^^^^^^^^^^^^^^

The shear strength is given by the normal stress :math:`\sigma_{yy}` multiplied by a constant friction coefficient :math:`\mu`

.. math::
   \tau^{str} = \sigma_{yy} \mu
   
Regularized Coulomb friction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The classical formulation of Coulomb friction might result in ill-posedness of the friction problem. A simplified regularization is thus implemented, based on `Prakash (1998) <https://asmedigitalcollection.asme.org/tribology/article/120/1/97/439195/Frictional-Response-of-Sliding-Interfaces>`_ and `Rubin and Ampuero (2007) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2006JB004337>`_ . The contact pressure is regularized as:

.. math::
   \frac{\tilde{\sigma}_{yy}}{dt} = -\frac{1}{t^*}\left(\tilde{\sigma}_{yy} - \sigma{yy} \right)

with :math:`t^*` a regularization parameter. The strength is computed as

.. math::
   \tau^{str} = \tilde{sigma}_{yy} \mu

Rate and state friction
-----------------------

The rate and state framework involve two functionnals: one for the friction coefficient :math:`f(v,\phi)` and one evolution law for the state variable :math:`\dot\phi = g(v,\phi)`. Several variations of the friction law and the evolution law are available in cRacklet. Their exact formulations are summarized here: 

Pure Weakening Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^

This formulation is based on the one originaly proposed by `Dieterich (1979) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB084iB05p02161>`_ and `Ruina (1983) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB088iB12p10359>`_ :

.. math::
   f(v,\phi) = f_0 + a \log \left(v/v_* \right) + b f_0 \log \left(\phi / \phi_* \right)

Standard Formulation
^^^^^^^^^^^^^^^^^^^^

This formulation is a modified version of the original formulation proposed by `Dieterich (1979) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB084iB05p02161>`_ and `Ruina (1983) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB088iB12p10359>`_, and has been proposed by `Bar-Sinai et al. (2012) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2011GL050554>`_.

.. math::
   f(v,\phi) = f_0 + a \log \left(1+v/v_* \right) + b f_0 \log \left(1 + \phi / \phi_* \right)

Regularized Formulation
^^^^^^^^^^^^^^^^^^^^^^^

This formulation is a generic N-shape friction law, as suggested by experimental observations (Bar-Sinai 2014). The steady-state friction is strenghtening at low and high velocities, and has a velocity-weakening branch at intermediate velocities, as observed for many materials in `Bar-Sinai et al. (2014) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2013JB010586>`_

.. math::
   f(v,\phi) = \left( 1 + b \log \left(1+ \frac{\phi}{\phi_*} \right) \right) \left( \frac{f_0}{  \sqrt{\left(  1+v_0^2 / v^2 \right)} } + a \log \left( 1+\frac{v}{v_*} \right) \right)

Regularized Weakening
^^^^^^^^^^^^^^^^^^^^^

This formulation is derived from the N-shape one, with the omission of the :math:`+1` term in the :math:`log` of :math:`\phi` term, resulting in the absence of the velocity-strengthening branch at high velocites.
   
.. math::
   f(v,\phi) = \left( 1 + b \log \left(\frac{\phi}{\phi_*} \right) \right) \left( \frac{f_0}{  \sqrt{\left(  1+v_0^2 / v^2 \right)} } + a \log \left( 1+\frac{v}{v_*} \right) \right)

Aging Law
^^^^^^^^^

Original evolution law proposed by `Dieterich (1979) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB084iB05p02161>`_.

.. math::
   g(\phi) = 1 - \frac{v \phi}{D}

Regularized Aging Law
^^^^^^^^^^^^^^^^^^^^^

Regularization of the original aging law.

.. math::
   g(\phi) = 1 - \frac{v \phi}{D} \sqrt{1 + \left(\frac{v_*}{v}\right)^2}

Slip Law
^^^^^^^^

Original slip law  proposed by `Ruina (1983) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB088iB12p10359>`_.

.. math::
   g(\phi) = - \frac{v \phi}{D} \log \left( \frac{v \phi}{D} \right)


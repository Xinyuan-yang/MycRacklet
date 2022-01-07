Interface Law
=============

cRacklet has several interfacial behavior implemented which are described in this section. These law are listed in :cpp:enum:`FractureLawType <FractureLawType>`.

Mixed-mode cohesive law
-----------------------

This cohesive law has its strength decreasing linearly with the opening. The opening considered here is the norm of the opening displacement, thus taking into account both normal and shear componenets. The critical displacement :math:`\delta_c` is the critical displacement required to transition from the peak stress :math:`\tau_c` to the residual value :math:`\tau_r` (By default is 0). The shear and normal components can be prescribed independantly. While the norm of the opening :math:`||\delta||` is lower than :math:`\delta_c`, the strength is given by:  

.. math::
   \tau^{str} = \left( \tau_c - \tau_r \right) \left(1- ||\delta||/\delta_c \right)

If :math:`||\delta||` is larger than :math:`\delta_c`, the strength is given by:

.. math::
   \tau^{str} = \tau_r

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
   \tau^{str} = \tilde{\sigma}_{yy} \mu

Rate and state friction
-----------------------

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
   
Custom interface law
--------------------

If you want to use a constitutive law that is not implemented in cRacklet, you will need to code it in C++. Please consider creating a merge request for your contribution in GitLab with appropriate documentation if you think that the constitutive behavior you implemented could be used by the community.


If you want to add an evolution law or a friction law that is part of the rate and state framework (depends on a state variable and the sliding velocity), you can define a simple *Struct* as done in the file *rate_and_state_formulations.hh*.

Adding a rate and state evolution law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For an evolution law, you need to define the evolution law :math:`g(\phi)` as the operator () and the evaluation of :math:`g(\phi)` in steady state as the method getSteadyState. In addition, you will need to create an initialization function for the class :cpp:class:`RateAndStateLaw <RateAndStateLaw>` that allocates a shared pointer to an object of your new "struct" (should be similar to :cpp:func:`initStateEvolution <RateAndStateLaw::initStateEvolution>`).

Adding a rate and state friction law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a friction law, you need to define the friction law :math:`f(v,\phi)` as the operator (), the derivative of the friction coefficient with respect to the velocity :math:`\partial f(v,\phi) / \partial v` as the method getTangent, the derivative of the steady-state friction coefficient with respect to the steady velocity :math:`\partial f(v_{ss}) / \partial v_{ss}` as the method getSteadyTangent, and a method that computes the state value satisfying steady state as getStableState.
In addition, an initialization function for the class :cpp:class:`RateAndStateLaw <RateAndStateLaw>` that allocates a shared pointer to an object of your new "struc" should be defined, similarly to :cpp:func:`initStandardFormulation <RateAndStateLaw::initStandardFormulation>`.

To create the interface itself, you should create a new entry in the enum :cpp:enum:`FractureLawType <FractureLawType>` with your new law, and create a templated version of the method :cpp:func:`createUniformInterface <Interfacer::createUniformInterface>` and :cpp:func:`createHeterogeneousInterface <Interfacer::createHeterogeneousInterface>` for the class :cpp:class:`Interfacer`. 

Creating a custom constitutive law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to define a constitutive law that is neither a cohesive law (strength depends only on opening) nor is defined in the rate and state framework (friction depends on velocity and a state variable), you will need to create a new class. This class should inherit from :cpp:class:`InterfaceLaw` and you need to define the following method:

- :cpp:func:`initInterfaceConditions <InterfaceLaw::initInterfaceConditions>`: computes the interface fields for the initial conditions.
- :cpp:func:`updateInterfaceConditions <InterfaceLaw::updateInterfaceConditions>`: computes the interface fields after the displacement has been updated according to the explicit time stepping and the dynamic stresses computed based on the displacement history.
- :cpp:func:`restart <InterfaceLaw::restart>`: save or load the fields required to restart a simulation. This should only concern additional fields used for the constitutive behavior implemented, as the displacement and velocity fields are handled by the mother class.

To solve for the interface fields, you can take inspiration from the classes :cpp:class:`CohesiveLaw` or :cpp:class:`RateAndStateLaw`.
With the class :cpp:class:`CohesiveLaw`, the system is well defined and the velocity is directly expressed in terms of the loading, the dynamic stresses and the material parameters. The normal and shear velocities are computed independently. The strength of the interface has to be evaluated against the stress to assess if interface opening is occuring.
The :cpp:class:`RateAndStateLaw` however uses an iterative Newton-Raphson procedure to find the value of the velocity and the state variable to satisfy that the interfacial stress is always equal to the strength.

To create the interface itself, you should create a new entry in the enum :cpp:enum:`FractureLawType <FractureLawType>` with your new law, and create a templated version of the method :cpp:func:`createUniformInterface <Interfacer::createUniformInterface>` and :cpp:func:`createHeterogeneousInterface <Interfacer::createHeterogeneousInterface>` for the class :cpp:class:`Interfacer`.

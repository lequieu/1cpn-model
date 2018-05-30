
1CPN in Lammps 
==================

Linking 1CPN with Lammps
-------------------------
Add Detailed instructions here


1CPN High Performance
-----------------------

mention the `fix balance` command

Custom 1CPN Potentials
-----------------------

While developing the 1CPN model, many custom potentials were necessary


Pair Zewdie
^^^^^^^^^^^^^^^^^^
Give functional form of potential?
And an explanation of the pair_coeffs arguments

*Usage:* ::

  pair_style zewdie ${pe000} ${pecc2} ${pe220} ${pe222} ${pe224} ${ps000} ${pscc2} ${ps220} ${ps222} ${ps224} 

  pair_coeff   1 1 zewdie ${pe0} ${ps0} ${cutoff}


* `pe000` - See :math:`\epsilon_{000}` in :eq:`zewdie-eps`
* `pecc2` - See :math:`\epsilon_{cc2}` in :eq:`zewdie-eps`
* `pe220` - See :math:`\epsilon_{220}` in :eq:`zewdie-eps`
* `pe222` - See :math:`\epsilon_{222}` in :eq:`zewdie-eps`
* `pe224` - See :math:`\epsilon_{224}` in :eq:`zewdie-eps`
* `ps000` - See :math:`\sigma_{000}` in :eq:`zewdie-sig`
* `pscc2` - See :math:`\sigma_{cc2}` in :eq:`zewdie-sig`
* `ps220` - See :math:`\sigma_{220}` in :eq:`zewdie-sig`
* `ps222` - See :math:`\sigma_{222}` in :eq:`zewdie-sig`
* `ps224` - See :math:`\sigma_{224}` in :eq:`zewdie-sig`
* `pe0` - See :math:`\epsilon_{0}` in :eq:`zewdie-eps`
* `ps0` - See :math:`\sigma_{0}` in :eq:`zewdie-sig`

.. math:: {U}_{Zewdie} ({r}_{ij},{\hat{f}}_i , {\hat{f}}_j \,;\, \sigma_0, \epsilon_0) = 
    4 \epsilon \left[ \left(\frac{\sigma_0}{r_{ij}-\sigma+\sigma_0}\right)^{12}- \left(\frac{\sigma_0}{r_{ij}-\sigma+\sigma_0)}\right)^6 \right]
    :label: zewdie

where 

.. math:: \sigma = \sigma_{0} [\sigma_{000}S_{000} + \sigma_{cc2}(S_{022} + S_{202}) + \sigma_{220}S_{220}  + \sigma_{222}S_{222}  + \sigma_{224}S_{224}]
   :label: zewdie-sig

.. math::   \epsilon = \epsilon_{0} [\epsilon_{000}S_{000} + \epsilon_{cc2}(S_{022} + S_{202}) + \epsilon_{220}S_{220} + \epsilon_{222}S_{222} + \epsilon_{224}S_{224}].
   :label: zewdie-eps

and the S-functions are defined as

.. math::

  S_{000} &= 1,\\
  S_{202} &= (3 a_1^2 - 1) / 2\sqrt{5},\\
  S_{022} &= (3 a_2^2 - 1) / 2\sqrt{5},\\
  S_{220} &= (3 a_0^2 - 1) / 2\sqrt{5},\\
  S_{222} &= \dfrac{1}{\sqrt{70}}(2 - 3 a_1^2 - 3 a_2^2 - 3 a_0^2 + 9 a_0 a_1 a_2),\\
  S_{224} &= \dfrac{1}{4\sqrt{70}}(1 + 2 a_0^2 - 5 a_1^2 - 5 a_2^2 - 20 a_0 a_1 a_2 + 35 a_1^2 a_2^2)\\
 
with

.. math::  a_0 &= {\hat{f}}_i \cdot {\hat{f}}_j, \\
           a_1 &= {\hat{f}}_i \cdot {\hat r}_{ij}, \\
           a_2 &= {\hat{f}}_j \cdot {\hat{r}}_{ij} \\


Discuss how sphere-ellipse interactions are handled



Pair Gauss Aniso
^^^^^^^^^^^^^^^^^^

*Usage:* ::

  pair_style gauss/aniso ${gauss_rcut}

  pair_coeff 2 3  sigma d0 r0 theta0 phi0 Ktheta Kphi

* `sigma` - Gaussian width. Given by :math:`\sigma` in :eq:`gauss`
* `d0` - Gaussian depth. Given by :math:`d_0` in :eq:`gauss`
* `r0` - Gaussian Center Position. Given by :math:`r_0` in :eq:`gauss`
* `theta0` - Equilibrium position with with respect to :math:`\theta`, where :math:`\theta = \arccos(\hat{r}_{ij} \cdot \hat{u}_i)`. See :math:`\theta_0` in :eq:`gauss-modulate`
* `phi0` - Equilibrium position with with respect to :math:`\phi`, where :math:`\phi = \arccos(\hat{r}_{ij} \cdot \hat{f}_i)`. See :math:`\phi_0` in :eq:`gauss-modulate`
* `Ktheta` - Width of modulating function with respect to :math:`\theta`. See :math:`K_\theta` in :eq:`gauss-modulate`
* `Kphi` - Width of modulating function with respect to :math:`\phi`. See :math:`K_\phi` in :eq:`gauss-modulate`

.. warning::
   In the implementation, of this potential, the atom_types of sites i and j cannot be the same.

   This is because the ith particle is always chosen to be the atom with the lower type index. 

   Not symmetric.

   For example 

.. math:: U_{gauss,aniso} &= f(K_\theta,\Delta \theta) f(K_\phi, \Delta \phi) \, {U}_{gauss} \nonumber \\
                  &= f(K_\theta,\Delta \theta) f(K_\phi, \Delta \phi) \left(-d_0 e^{-(r-r_0)^2 / 2\sigma^2}\right)
  :label: gauss



.. math:: f(K_\theta,\Delta \theta) = 
  \begin{cases}
      1                                       &-\frac{\pi}{2K_\theta} < \Delta \theta < \frac{\pi}{2K_\theta} \\
      1-\cos^2\left( K_\theta \Delta \theta \right)  & \frac{-\pi}{K_\theta}<\Delta \theta < \frac{-\pi}{2K_\theta} \textrm{ or } \frac{\pi}{2K_\theta} < \Delta \theta < \frac{\pi}{K_\theta}\\
      0                                       &\Delta \theta < -\frac{\pi}{K_\theta} \textrm{ or } \Delta \theta > \frac{\pi}{K_\theta} 
  \end{cases}
  :label: gauss-modulate


Angle Orient
^^^^^^^^^^^^^^^^^^

*Usage:*:: 

  angle_style orient
  angle_coeff 1 angle_{f,v,u} ktheta1 ktheta2 kphi theta1 theta2 phi

* `angle_vector` - Possible values `angle_f`, `angle_v`, or `angle_u`. Defines whether :math:`\hat{w} = \{\hat{f},\hat{v},\hat{u}\}`. 
* `ktheta1` - Spring constant for deformations in :math:`\theta_1`. See :eq:`angle-orient`.
* `ktheta2` - Spring constant for deformations in :math:`\theta_2`. See :eq:`angle-orient`.
* `kphi` - Spring constant for deformations in :math:`\phi`. See :eq:`angle-orient`.
* `theta1` - Equilibrium value of :math:`\theta_1`. See :math:`\theta_{1,0}` in :eq:`angle-orient`.
* `theta2` - Equilibrium value of :math:`\theta_2`. See :math:`\theta_{2,0}` in :eq:`angle-orient`.
* `phi` - Equilibrium value of :math:`\phi`. See :math:`\phi_{0}` in :eq:`angle-orient`.


.. math:: U = \frac{1}{2} \left( k_{\theta_1} (\theta_1 - \theta_{1,0})^2 + 
    k_{\theta_2} (\theta_2 - \theta_{2,0})^2 + 
    k_{\phi} (\phi - \phi_{0})^2 \right)
    :label: angle-orient

where :math:`\theta_1, \theta_2, \phi` are given by

.. math:: \theta_1 &= \arccos (\hat{w}_i \cdot \hat{r}_{ij}) \\
          \theta_2 &= \arccos (\hat{w}_j \cdot \hat{r}_{ij}) \\
          \phi &= \arccos (\hat{w}_i \cdot \hat{w}_j)
 



Angle Orient Cosine
^^^^^^^^^^^^^^^^^^^^

*Usage:*:: 

  angle_style orient/cosine
  angle_coeff 1 angle_{f,v,u} ktheta1 ktheta2 kphi theta1 theta2 phi

* `angle_vector` - Possible values `angle_f`, `angle_v`, or `angle_u`. Defines whether :math:`\hat{w} = \{\hat{f},\hat{v},\hat{u}\}`. 
* `ktheta1` - Spring constant for deformations in :math:`\theta_1`. See :eq:`angle-orient-cos`.
* `ktheta2` - Spring constant for deformations in :math:`\theta_2`. See :eq:`angle-orient-cos`.
* `kphi` - Spring constant for deformations in :math:`\phi`. See :eq:`angle-orient-cos`.
* `theta1` - Equilibrium value of :math:`\theta_1`. See :math:`\theta_{1,0}` in :eq:`angle-orient-cos`.
* `theta2` - Equilibrium value of :math:`\theta_2`. See :math:`\theta_{2,0}` in :eq:`angle-orient-cos`.
* `phi` - Equilibrium value of :math:`\phi`. See :math:`\phi_{0}` in :eq:`angle-orient-cos`.


.. math:: U = \frac{1}{2} \left[ k_{\theta_1} \left(1 - \cos(\theta_1 - \theta_{1,0})\right) + 
          k_{\theta_2} \left(1 - \cos(\theta_2 - \theta_{2,0})\right) + 
          k_{\phi} \left(1 - \cos(\phi - \phi_{0})\right) \right]
    :label: angle-orient-cos

where :math:`\theta_1, \theta_2, \phi` are given by

.. math:: \theta_1 &= \arccos (\hat{w}_i \cdot \hat{r}_{ij}) \\
          \theta_2 &= \arccos (\hat{w}_j \cdot \hat{r}_{ij}) \\
          \phi &= \arccos (\hat{w}_i \cdot \hat{w}_j)
 


Angle WLC Twist
^^^^^^^^^^^^^^^^^^

*Usage:*:: 

  angle_style wlctwist
  angle_coeff 2 wlctwist ${kalign} ${ktwist} ${omega0} 

* `kalign` - Alignment spring constant
* `ktwist` - Twist spring constant
* `omega0` - Equilibrium Twist 

This potential is based off of the implementation of Brackley et al.

.. math:: U = k_{\omega} \left(1-\cos(\omega_i - \omega_{0})\right) + k_{\psi} \left(1-\cos (\psi_i)\right)
    :label: wlctwist


Maintaing Compatability with LAMMPS
------------------------------------

In order to maintain compatability of 1CPN with the most recent version of LAMMPS it is helpful to know 
which core LAMMPS potential the 1CPN potentials were derived from. By seeing what changed between the core 
LAMMPS potentials (i.e. diff old and new version), it is typically straightforward to make the necessary 
minor changes to the 1CPN potential to allow 1CPN-LAMMPS to compile.

When trying a new version of LAMMPS, be sure to run the integration tests in `${D_1CPN}/test/integ_tests`, to make sure the model is behaving correctly

  ============================  ================================
        1CPN Potential             Original Lammps Potential
  ============================  ================================
  pair_style zewdie             pair_style gayberne
  pair_style gauss/aniso        pair_style gayberne, pair_style gauss
  angle_style wlctwist          angle_style wlctwist (Brackley2014)
  angle_style orient            angle_style wlctwist
  angle_style orient/cosine     angle_style wlctwist
  ============================  ================================



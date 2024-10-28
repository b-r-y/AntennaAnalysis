"""Pyramidal horn class."""
import sys
import os
import numpy as np
import numpy.typing as npt
import scipy.constants as spc
from scipy.optimize import newton
from scipy.special import fresnel

from .common import CommonPattern
from ..includes.utilities import get_p_struct, PatternType


class PyramidalHorn(CommonPattern):
    """Analytical pyramidal horn class.

    Calculates approximate and exact solutions to the pyramidal horn and
    generates analytical radiation pattern.

    All calculations and definitions of variables according to [BALANIS]_.
    See :numref:`fig_pyramidal_horn`

    .. _fig_pyramidal_horn:
    .. figure:: ../img/pyramidal_horn.png
        :align: center
        :width: 600px

        Pyramidal horn definitions. Figure copied from [BALANIS]_ for
        completeness.

    Args:
        frequency (float): frequency of the horn in [Hz]
        a (float): larger side of the rectangular feed waveguide in [m]
        b (float): smaller side of the rectangular feed waveguide in [m]
        a1 (float): larger side of the horn open aperture in [m]
        b1 (float): smaller side of the horn open aperture in [m]
        rho1 (float): E plane radius as shown in \
            :numref:`fig_pyramidal_horn` in [m]
        rho2 (float): H plane radius as shown in \
            :numref:`fig_pyramidal_horn` in [m]
        resolution (:obj:`list` of two :obj:`int`): theta and phi steps for\
            pattern resolutions in degrees
        efficiency (float): total pattern efficiency. The pattern will be \
            normalized to this value. It is defined in percept from 0 to 100.
    """

    def __init__(self,
                 frequency: float,
                 a: float, b: float,
                 a1: float, b1: float,
                 rho1: float, rho2: float,
                 resolution: list[int],
                 efficiency: float):
        """Initilize pyramidal horn class."""
        self.resolution = resolution
        self.efficiency = efficiency
        self.frequency = frequency
        self._a = a
        self._b = b
        self._a1 = a1
        self._b1 = b1
        self._rho1 = rho1
        self._rho2 = rho2
        self._pattern: PatternType = get_p_struct()

    def export_mechanical_design_parameters(self) -> str:
        """Export mechanical parameters of the horn.

        Calculates all mechanical properties of the horn and:

            * exports a Solidworks formated file for parametrized structure \
                in the output folder
            * returns a string that simply lists those mechanical values in \
                human readable form

        Returns:
            str: string with human readable mechanical dimensions according to\
                :numref:`fig_pyramidal_horn`
        """
        # Exporting SW design file
        if not os.path.exists('output'):
            os.makedirs('output')
        eq_f = open(os.path.join('output', "equations.txt"), "w",
                    encoding="UTF-8")
        eq_f.write((f'"rho_1" = {self.rho1*1000:.2f}\n'
                    f'"b1" = {self.b1*1000:.2f}\n'
                    f'"a1" = {self.a1*1000:.2f}\n'
                    f'"b" = {self.b*1000:.2f}\n'
                    f'"a" = {self.a*1000:.2f}\n'
                    f'"wall_thickness" = 2\n'
                    f'"feed_length" = 10\n'
                    f'"rho_e" = sqr ( ( ( "b1" / 2 ) ^ 2 ) + "rho_1" ^ 2 )\n'
                    f'"PE"= ( "b1" - "b" ) * '
                    f'sqr ( ( ( "rho_e" / "b1" ) ^ 2 ) - 0.25 )'
                    )
                   )
        eq_f.close()
        string = (
                  f'"a" = {self.a*1000:.2f} [mm]\n'
                  f'"b" = {self.b*1000:.2f} [mm]\n'
                  f'"b1" = {self.b1*1000:.2f} [mm]\n'
                  f'"a1" = {self.a1*1000:.2f} [mm]\n'
                  f'"rho1" = {self.rho1*1000:.2f} [mm]\n'
                  f'"rho2" = {self.rho2*1000:.2f} [mm]\n'
                  f'"rhoE" = {self.rhoE*1000:.2f} [mm]\n'
                  f'"rhoH" = {self.rhoH*1000:.2f} [mm]\n'
                  f'"pE" = {self.pE*1000:.2f} [mm]\n'
                  f'"pH" = {self.pH*1000:.2f} [mm]\n'
                  f'"psiE" = {self.psiE*(180/np.pi):.2f} [deg]\n'
                  f'"psiH" = {self.psiH*(180/np.pi):.2f} [deg]\n'
                  )
        return string

    def calculate_approximate_directivity(self) \
            -> tuple[float, float, float]:  # noqa: D301
        """Calculate directivities for E, H pyramidal horn solutions.

        Calculates approximate directivities for E, H and pyramidal horn
            solutions using equations 13-19, 13-41 and 13-52a respectively
            from [BALANIS]_.

        .. math::
            D_{E} = \\frac{64a\\rho_{1}}{\\pi\\lambda b_{1}} \
                \\left[ \
                C^2 \
                \\left( \\frac{b_{1}}{\\sqrt{2\\lambda\\rho_{1}}} \\right) \
                + S^2 \
                \\left( \\frac{b_{1}}{\\sqrt{2\\lambda\\rho_{1}}} \\right) \
                \\right]

        .. math::
            D_{H} = \\frac{4b\\pi\\rho_{2}}{a_{1}\\lambda} \
                \\times \
                \\left\\{ \
                \\left[ C(u) - C(v) \\right ]^2 + \
                \\left[ S(u) - S(v) \\right ]^2 \
                \\right\\}

        where

        .. math::
            u &= \\frac{1}{\\sqrt{2}} \
                \\left( \
                \\frac{\\sqrt{\\lambda\\rho_{2}}}{a_{1}} + \
                \\frac{a_{1}}{\\sqrt{\\lambda\\rho_{2}}} \
                \\right)

            v &= \\frac{1}{\\sqrt{2}} \
                \\left( \
                \\frac{\\sqrt{\\lambda\\rho_{2}}}{a_{1}} - \
                \\frac{a_{1}}{\\sqrt{\\lambda\\rho_{2}}} \
                \\right)

        and :math:`C(x)` and :math:`S(x)` are the Fresnel integrals.

        .. math::
            D_{P} = \\frac{\\pi\\lambda^2}{32ab}D_{E}D_{H}

        Returns:
            tuple: the values of DE, DH and DP in linear format
        """
        # E plane
        fresnel_argument_e_plane = self.b1 / np.sqrt(2 *
                                                     self.rho1 *
                                                     self.wavelength0)
        S, C = fresnel(fresnel_argument_e_plane)
        D_E = ((64*self.a*self.rho1) / (np.pi*self.b1*self.wavelength0)) * \
              ((C**2) + (S**2))
        # H plane
        u = 1/np.sqrt(2) * \
            ((np.sqrt(self.rho2*self.wavelength0)/self.a1) +
             (self.a1/np.sqrt(self.rho2*self.wavelength0)))
        v = 1/np.sqrt(2) * \
            ((np.sqrt(self.rho2*self.wavelength0)/self.a1) -
             (self.a1/np.sqrt(self.rho2*self.wavelength0)))
        S_u, C_u = fresnel(u)
        S_v, C_v = fresnel(v)
        D_H = ((4*np.pi*self.b*self.rho2)/(self.a1*self.wavelength0)) * \
              (((C_u - C_v)**2) + ((S_u-S_v)**2))
        D_P = ((np.pi*self.wavelength0**2)/(32*self.a*self.b))*D_E*D_H
        return D_E, D_H, D_P

    def get_pattern(self) -> PatternType:
        """Generate analytical radiation pattern from mechanical design values.

        Generation is according to section 13.4 in [BALANIS]_.

        Returns:
            dict: pattern according to the structure \
                :func:`~antenna_analysis.includes.utilities.get_p_struct`
        """
        vth_mesh, vph_mesh = self._get_mesh()

        # E plane
        ky = self.wavenumber0*np.sin(vth_mesh)*np.sin(vph_mesh)
        t1 = np.sqrt(1/(np.pi*self.wavenumber0*self.rho1)) * \
            (-((self.wavenumber0*self.b1)/2) - (ky*self.rho1))
        t2 = np.sqrt(1/(np.pi*self.wavenumber0*self.rho1)) * \
            (((self.wavenumber0*self.b1)/2) - (ky*self.rho1))

        S_t1, C_t1 = fresnel(t1)
        S_t2, C_t2 = fresnel(t2)
        I2 = np.sqrt((np.pi*self.rho1)/self.wavenumber0) * \
            np.exp(1j*((self.rho1*ky**2)/(2*self.wavenumber0))) * \
            ((C_t2 - C_t1) - 1j*(S_t2-S_t1))
        # H plane
        kx_prime = self.wavenumber0 * np.sin(vth_mesh) * np.cos(vph_mesh) + \
            np.pi/self.a1
        t1_prime = np.sqrt(1/(np.pi*self.wavenumber0*self.rho2)) * \
            (-((self.wavenumber0*self.a1)/2)-kx_prime*self.rho2)
        t2_prime = np.sqrt(1/(np.pi*self.wavenumber0*self.rho2)) * \
            (((self.wavenumber0*self.a1)/2)-kx_prime*self.rho2)

        kx_second = self.wavenumber0 * np.sin(vth_mesh) * np.cos(vph_mesh) - \
            np.pi/self.a1
        t1_second = np.sqrt(1/(np.pi*self.wavenumber0*self.rho2)) * \
            (-((self.wavenumber0*self.a1)/2)-kx_second*self.rho2)
        t2_second = np.sqrt(1/(np.pi*self.wavenumber0*self.rho2)) * \
            (((self.wavenumber0*self.a1)/2)-kx_second*self.rho2)

        S_t1_prime, C_t1_prime = fresnel(t1_prime)
        S_t2_prime, C_t2_prime = fresnel(t2_prime)
        S_t1_second, C_t1_second = fresnel(t1_second)
        S_t2_second, C_t2_second = fresnel(t2_second)
        I1 = 0.5 * np.sqrt((np.pi * self.rho2) / self.wavenumber0) * \
            (np.exp(1j*((self.rho2*kx_prime**2)/(2*self.wavenumber0))) *
             ((C_t2_prime - C_t1_prime) - 1j*(S_t2_prime - S_t1_prime))
             +
             np.exp(1j*((self.rho2*kx_second**2)/(2*self.wavenumber0))) *
             ((C_t2_second-C_t1_second) - 1j*(S_t2_second-S_t1_second))
             )

        e_th: npt.NDArray[np.complex64] = 1j * \
            ((self.wavenumber0 * np.exp(-1j * self.wavenumber0)) /
             (4 * np.pi)) * \
            (np.sin(vph_mesh) * (1 + np.cos(vth_mesh)) * I1 * I2)
        e_ph: npt.NDArray[np.complex64] = 1j * \
            ((self.wavenumber0 * np.exp(-1j * self.wavenumber0)) /
             (4 * np.pi)) * \
            (np.cos(vph_mesh) * (np.cos(vth_mesh) + 1) * I1 * I2)

        e_th = np.transpose(e_th)
        e_ph = np.transpose(e_ph)

        e_th, e_ph = self._normalize(e_th, e_ph)

        # Store the generated pattern
        self._pattern['source'] = 'Analytical Pyramidal Horn'
        self._pattern['file'] = 'N/A'
        self._pattern['frequency'] = self.frequency
        self._pattern['port'] = 1
        self._pattern['th_step'] = self._resolution[0]
        self._pattern['ph_step'] = self._resolution[1]
        self._pattern['ff'] = np.stack((e_th, e_ph))

        return self._pattern

    @property
    def a(self) -> float:
        return self._a

    @property
    def b(self) -> float:
        return self._b

    @property
    def a1(self) -> float:
        return self._a1

    @property
    def b1(self) -> float:
        return self._b1

    @property
    def rho1(self) -> float:
        return self._rho1

    @property
    def rho2(self) -> float:
        return self._rho2

    @property
    def rhoE(self) -> float:
        return np.sqrt(((self.b1/2)**2) + self.rho1**2)

    @property
    def rhoH(self) -> float:
        return np.sqrt(((self.a1/2)**2) + self.rho2**2)

    @property
    def pE(self) -> float:
        return (self.b1 - self.b)*np.sqrt(((self.rhoE/self.b1)**2) - (1/4))

    @property
    def pH(self) -> float:
        return (self.a1 - self.a)*np.sqrt(((self.rhoH/self.a1)**2) - (1/4))

    @property
    def psiE(self) -> float:
        return np.arcsin(self.b1/(2*self.rhoE))

    @property
    def psiH(self) -> float:
        return np.arcsin(self.a1/(2*self.rhoH))


class PyramidalHornBuilder:
    """Build a pyramidal horn object."""

    def __init__(self):
        """Initialize the builder object."""
        self._instance = None

    def __call__(self, frequency: float,
                 a: float, b: float,
                 target_gain: float = 15,
                 a1: float = 0, b1: float = 0,
                 rho1: float = 0, rho2: float = 0,
                 resolution: list[int] = None,
                 efficiency: float = 100,
                 **_ignored) -> PyramidalHorn:
        """Build a pyramidal horn object.

        The function first checks if the given `frequency` is supported by
        the dimensions of the waveguide given by `a` and `b`.
        If not it terminates the program.

        If any of `a1`, `b1`, `rho1` or `rho2` are left as default
        `None`, `target_gain` is used
        to calculate the optimum horn according to 13.4.3 in [BALANIS]_.
        Otherwise the submitted
        values are used and `target_gain` is ignored.

        Args:
            frequency (float): frequency of the horn in Hz
            a (float): larger side of the rectangular feed waveguide in [m]
            b (float): smaller side of the rectangular feed waveguide in [m]
            target_gain (float): target design gain for cases where not all\
                mechanical properties are given. Value is in dB.
            a1 (float): larger side of the horn open aperture in [m]
            b1 (float): smaller side of the horn open aperture in [m]
            rho1 (float): E plane radius as shown in \
                :numref:`fig_pyramidal_horn` in [m]
            rho2 (float): H plane radius as shown in \
                :numref:`fig_pyramidal_horn` in [m]
            resolution (:obj:`list` of two :obj:`int`): \
                theta and phi steps for pattern resolutions in degrees
            efficiency (float): total pattern efficiency. \
                The pattern will be normalized to this value. \
                It is defined in percent from 0 to 100.

        Returns:
            object: an object of class \
                :class:`~antenna_analysis.imports.pyramidal_horn.PyramidalHorn`

        """
        if resolution is None:
            resolution = [1, 1]
        # self._frequency = frequency
        wavelength = spc.speed_of_light/frequency
        # From Eq. 3.84 in Pozar
        te10_cutoff_frequency = 1/(2*a*np.sqrt(spc.mu_0*spc.epsilon_0))
        te20_cutoff_frequency = 2/(2*a*np.sqrt(spc.mu_0*spc.epsilon_0))
        if te10_cutoff_frequency < frequency < te20_cutoff_frequency:
            if any([x == 0 for x in [a1, b1, rho1, rho2]]):
                a1, b1, rho1, rho2 = self.__calculate_optimal_pyramidal_horn(
                        target_gain, a, b, wavelength)
            self._instance = PyramidalHorn(frequency, a, b,
                                           a1, b1, rho1, rho2,
                                           resolution, efficiency)
        else:
            print(f'Requested frequency NOT supported by waveguide!\n'
                  f'TE10 cutoff frequency: '
                  f'{te10_cutoff_frequency/1E6:.2f} MHz\n'
                  f'Requested frequency: {frequency/1E6:.2f} MHz\n'
                  f'TE20 cutoff frequency: '
                  f'{te20_cutoff_frequency/1E6:.2f} MHz\n'
                  f'Calculation based on '
                  f'a={a*1000} and b={b*1000} mm size waveguide\n'
                  f'Correct the input before continuing!\n'
                  f'Program will now exit...')
            sys.exit()
        return self._instance

    def __calculate_optimal_pyramidal_horn(self, target_gain_db: float,
                                           a: float, b: float,
                                           wavelength: float) \
            -> tuple[float, float, float, float]:
        target_gain_lin = 10**(target_gain_db/10)  # linear target gain
        # calculate TE10 mode cutoff frequency
        chi: float = target_gain_lin/(2*np.pi*np.sqrt(2*np.pi))
        chi = newton(
            self.__optimization_subroutine, chi,
            args=(a, b, wavelength, target_gain_lin, ),
            tol=np.finfo(float).eps)
        rho_e: float = chi * wavelength
        rho_h = ((target_gain_lin**2)*wavelength)/(8*chi*np.pi**3)
        a1 = np.sqrt(3*wavelength*rho_h)
        b1 = np.sqrt(2*wavelength*rho_e)
        p_e = (b1 - b)*np.sqrt(((rho_e/b1)**2) - (1/4))
        p_h = (a1 - a)*np.sqrt(((rho_h/a1)**2) - (1/4))
        psi_e = np.arcsin(b1/(2*rho_e))
        psi_h = np.arcsin(a1/(2*rho_h))
        assert round(p_e, 6) == round(p_h, 6)
        rho1 = rho_e*np.cos(psi_e)
        rho2 = rho_h*np.cos(psi_h)
        return a1, b1, rho1, rho2

    def __optimization_subroutine(self, x: float,
                                  a: float,
                                  b: float,
                                  wavelength: float,
                                  target_gain_lin: float) -> float:
        # implements Equation 13-56 from Balanis
        left_side = ((np.sqrt(2 * x) - b / (wavelength))**2) * (2 * x - 1)
        right_side = (((target_gain_lin / (2 * np.pi * np.sqrt(x))) *
                      (np.sqrt(3 / (2 * np.pi))) - a / (wavelength))**2) * \
                     (((target_gain_lin**2)/(6 * x * np.pi**3)) - 1)
        return left_side - right_side

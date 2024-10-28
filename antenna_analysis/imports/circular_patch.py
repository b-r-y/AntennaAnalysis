"""Analytical circular patch."""
import sys
from typing import Any
import numpy as np
import numpy.typing as npt
import scipy.constants as spc
import scipy.special as sps

from .common import CommonPattern
from ..includes.utilities import get_p_struct, PatternType


class CircularPatch(CommonPattern):
    """Analytical circular patch class.

    Generates analytical radiation pattern of a circular patch antenna.

    All calculations and definitions of variables according to [BALANIS]_.
    See :numref:`fig_circular_patch`

    .. _fig_circular_patch:
    .. figure:: ../img/circular_patch.png
        :align: center
        :width: 600px

        Circular patch definitions. Figure copied from [BALANIS]_ for
        completeness.
    """

    def __init__(self, frequency: float,
                 substrate_height: float,
                 epsilon_r: float,
                 patch_radius: float,
                 resolution: list[int],
                 efficiency: float):
        """Analytical circular patch parameters.

        Args:
            frequency (float): frequency of the antenna in [Hz]
            substrate_height (float): substrate height in [mm]. Internally \
                this is converted to and used in meters.
            epsilon_r (float): substrate dielectric constant - dimenssionless
            patch_radius (float): patch radius in [mm]. Internally this is\
                converted to and used in meters.
            resolution (:obj:`list` of two :obj:`int`): theta and phi steps\
                for pattern resolutions in degrees
            efficiency (float): total pattern efficiency. The pattern will be \
                normalized to this value. It is defined in percept from 0 to\
                100.
        """
        self.resolution = resolution
        self.efficiency = efficiency
        self.frequency = frequency
        self.substrate_height = substrate_height
        self.epsilon_r = epsilon_r
        self.patch_radius = patch_radius
        self._resolution = resolution
        self._pattern: PatternType = get_p_struct()

    def describe(self) -> str:
        """Return a description string.

        Generates and returns a description string of the antenna.

        Returns:
            str: string with human readable description of dimensions\
                according to :numref:`fig_circular_patch`
        """
        string = (
            f'Analytical circular patch antenna:\n'
            f'Frequency (f0) = {self.frequency/1E6:.2f} [MHz]\n'
            f'Substrate height (h) = {self.substrate_height:.2f} [mm]\n'
            f'Substrate dielectric constant (e_r) = {self.epsilon_r:.2f} [-]\n'
            f'Patch radius (a) = {self.patch_radius:.2f} [mm]\n'
            f'Efficiency = {self.efficiency:.2f} [%]\n'
        )
        return string

    def get_pattern(self) -> PatternType:
        """Generate analytical radiation pattern from mechanical values.

        Generation is according to section 14.3 in [BALANIS]_.

        Returns:
            dict: pattern according to the structure\
                :func:`~antenna_analysis.includes.utilities.get_p_struct`
        """
        a = self._patch_radius
        h = self._substrate_height
        e_r = self.epsilon_r
        V0 = 1
        r = 1

        vth, vph = self._get_vectors()
        vth_mesh, vph_mesh = self._get_mesh()

        # BALANIS Eq. 14-67
        a_e = a*np.sqrt(1 + ((2*h)/(np.pi*a*e_r)) *
                        (np.log((np.pi*a)/(2*h)) + 1.7726))

        bessel_arg = self.wavenumber0*a_e*np.sin(vth)
        J02P = sps.jv(0, bessel_arg) - sps.jv(2, bessel_arg)
        J02 = sps.jv(0, bessel_arg) + sps.jv(2, bessel_arg)
        constant_prefix = 1j * \
            (self.wavenumber0 * a_e * V0 * np.exp(-1j*self.wavenumber0*r)) / \
            (2*r)
        e_th: npt.NDArray[np.complex64] = \
            - constant_prefix * J02P * np.cos(vph_mesh)
        e_th = np.transpose(e_th)
        e_th[vth.__len__()//2 + 1:, :] = \
            np.zeros((vth.__len__()//2, vph.__len__()-1), dtype=complex)
        e_ph: npt.NDArray[np.complex64] = \
            constant_prefix * J02 * np.cos(vth_mesh) * np.sin(vph_mesh)
        e_ph = np.transpose(e_ph)
        e_ph[vth.__len__()//2 + 1:, :] = \
            np.zeros((vth.__len__() // 2, vph.__len__()-1), dtype=complex)

        e_th, e_ph = self._normalize(e_th, e_ph)

        # Store the generated pattern
        self._pattern['source'] = 'Analytical Circular Patch'
        self._pattern['file'] = 'N/A'
        self._pattern['frequency'] = self.frequency
        self._pattern['port'] = 1
        self._pattern['th_step'] = self._resolution[0]
        self._pattern['ph_step'] = self._resolution[1]
        self._pattern['ff'] = np.stack((e_th, e_ph))

        return self._pattern

    @property
    def substrate_height(self) -> float:
        """Substrate height in [mm]. Stored and used in meters."""
        return self._substrate_height*1000

    @substrate_height.setter
    def substrate_height(self, value: float) -> None:
        """Set the substrate height in [mm]. Stored and used in meters."""
        self._substrate_height = value/1000

    @property
    def epsilon_r(self) -> float:
        """Dielectric constant [dimensionless]."""
        return self._epsilon_r

    @epsilon_r.setter
    def epsilon_r(self, value: float) -> None:
        self._epsilon_r = value

    @property
    def patch_radius(self) -> float:
        """Patch radius in [mm]. Stored and used in meters."""
        return self._patch_radius*1000

    @patch_radius.setter
    def patch_radius(self, value: float) -> None:
        """Patch radius in [mm]. Stored and used in meters."""
        self._patch_radius = value/1000


class CircularPatchBuilder:
    """Builds a circular patch object."""

    def __init__(self):
        """Initialize Circular patch builder."""
        self._instance = None

    def __call__(self, frequency: float = 1E9,
                 substrate_height: float = 10,
                 epsilon_r: float = 1,
                 patch_radius: float = 0,
                 resolution: list[int] = None,
                 efficiency: float = 100,
                 **_ignored: Any) -> CircularPatch:  # noqa: D301
        """Build a circular patch object.

        The function first checks if the given `substrate_height` is
        significantly smaller than the wavelength at the `frequency` for the
        given `epsilon_r`:

        .. math::
            \\frac{h}{\\lambda} < 0.5

        If not it terminates the program.

        Then it processes the patch radius and if the default 0 is not
        overwritten, calculates the ideal patch radius according to [BALANIS]_
        Eq. 14-69

        Args:
            frequency (float): frequency of the antenna in [Hz]
            substrate_height (float): substrate height in [mm]. Internally
                this is converted to and used in meters.
            epsilon_r (float): substrate dielectric constant - dimenssionless
            patch_radius (float): patch radius in [mm]. Internally this is
                converted to and used in meters.
            resolution (:obj:`list` of two :obj:`int`): theta and phi steps
                for pattern resolutions in degrees
            efficiency (float): total pattern efficiency. The pattern will be
                normalized to this value. It is defined in percept from
                0 to 100.

        Returns:
            object: an object of class
                :class:`~antenna_analysis.imports.circular_patch.CircularPatch`
        """
        if resolution is None:
            resolution = [1, 1]
        substrate_heigth_to_wavelength_ratio: float = \
            (substrate_height/1000) / \
            (spc.speed_of_light/(frequency*np.sqrt(epsilon_r)))
        ratio = 10
        # wavelength is twice as big
        if substrate_heigth_to_wavelength_ratio > 1/ratio:
            print(
                f'Substrate height ({substrate_height/1000:.4f} [m])'
                f' to wavelength ('
                f'{spc.speed_of_light/(frequency*np.sqrt(epsilon_r)):.4f} [m])'
                f' ratio not met!\n'
                f'Substrate heigh should be at least {ratio} times smaller '
                f'than the wavelength but {ratio}x{substrate_height/1000:.4f}='
                f'{ratio*substrate_height/1000:.4f} > '
                f'{spc.speed_of_light/(frequency*np.sqrt(epsilon_r)):.4f}\n'
                f'Correct the input before continuing!\n'
                f'Program will now exit...'
            )
            sys.exit()

        if patch_radius <= 0:
            # if the patch radius is not specified it is automatically
            # calculated from the frequency this calculation is done in
            # centimeters and the result is converted to meters
            F = (8.791*1E9) / (frequency * np.sqrt(epsilon_r))
            a = F / np.sqrt(1 + 2 * (substrate_height / 10) /
                            (np.pi * epsilon_r * F) *
                            (
                              np.log(np.pi * F / (2 * (substrate_height / 10)))
                              + 1.7726
                            )
                            )
            patch_radius = a*10

        self._instance = CircularPatch(frequency, substrate_height, epsilon_r,
                                       patch_radius, resolution, efficiency)
        return self._instance

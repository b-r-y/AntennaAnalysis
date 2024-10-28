"""Common antenna pattern class."""

from abc import ABC
import numpy as np
import numpy.typing as npt
import scipy.constants as spc


class CommonPattern(ABC):
    """Defines basic properties of antenna patterns."""

    @property
    def wavelength0(self) -> float:  # noqa: D301
        """Wavelength in vacuum.

        Wavelength in vacuum calculated as:

            :math:`\\lambda_{0} = \\frac{speed\\, of\\, light}{frequency}`

        Speed of light value is for vacuum from the scipy constants library.

        Returns:
            float: wavelength in vacuum in :math:`[m]`
        """
        return spc.speed_of_light/self.frequency

    @property
    def wavenumber0(self) -> float:  # noqa: D301
        """Wave number in vacuum.

        Wave number in vacuum calculated as:

            :math:`k_{0} = \\frac{2\\pi}{\\lambda_{0}}`

        See :meth:`~imports.common.CommonPattern.wavelength0` for details.

        Returns:
            float: wave number in vacuum in :math:`[m^-1]`
        """
        return (2*np.pi)/self.wavelength0

    @property
    def resolution(self) -> list[int]:
        """Pattern theta and phi resolution.

        Returns:
            :obj:`list` of two :obj:`int`: theta and phi steps\
                for pattern resolutions in degrees
        """
        return self._resolution

    @resolution.setter
    def resolution(self, resolution: list[int]) -> None:
        self._resolution = resolution

    @property
    def efficiency(self) -> float:  # noqa D301
        """Pattern efficiency in percent.

        Returns:
            float: efficiency in :math:`\\%`
        """
        return self._efficiency

    @efficiency.setter
    def efficiency(self, efficiency: float) -> None:
        self._efficiency = efficiency

    @property
    def frequency(self) -> float:
        """Pattern frequency in Hertz.

        Returns:
            float:frequency in :math:`[Hz]`
        """
        return self._frequency

    @frequency.setter
    def frequency(self, frequency: float) -> None:
        self._frequency = frequency

    def _get_vectors(self) -> tuple[npt.NDArray[np.float64],
                                    npt.NDArray[np.float64]]:
        vth = np.array(range(0, 180+self.resolution[0], self.resolution[0]))\
                * np.pi/180
        # vph is defined 0:360 for the integration only. For normal matrix
        # storage the 360 degrees is not needed
        vph = np.array(range(0, 360+self.resolution[1], self.resolution[1]))\
            * np.pi/180
        return vth, vph

    def _get_mesh(self) -> list[npt.NDArray[np.float64]]:
        vth, vph = self._get_vectors()
        return np.meshgrid(vth, vph[0:-1])

    def _normalize(self,
                   e_th: npt.NDArray[np.complex64],
                   e_ph: npt.NDArray[np.complex64]) -> \
            tuple[npt.NDArray[np.complex64],
                  npt.NDArray[np.complex64]]:
        vth, vph = self._get_vectors()
        # Normalization
        # compute current total power per point
        g_abs = np.abs(e_th)**2 + np.abs(e_ph)**2
        for i in range(0, int(180//self.resolution[0]+1)):
            g_abs[i, :] = g_abs[i, :] * np.sin(vth[i])
        # Close the sphere:
        g_abs = np.append(g_abs, np.transpose([g_abs[:, 0]]), axis=1)
        g_norm = (self.efficiency/100) * \
                 (
                 (4*np.pi)/(np.trapz(np.trapz(g_abs, vth, axis=0), vph))
                 )
        e_th = e_th*np.sqrt(g_norm)
        e_ph = e_ph*np.sqrt(g_norm)

        return e_th, e_ph

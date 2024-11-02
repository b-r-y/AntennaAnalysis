"""3GPP antenna element."""

import sys
from dataclasses import dataclass, field
import numpy as np
import numpy.typing as npt
from antenna_analysis.imports.common import CommonPattern
from antenna_analysis.includes.utilities import get_p_struct, PatternType


@dataclass
class TGPPAntenna(CommonPattern):
    """3GPP antenna models class.

    This class generates antennas according to 3GPP.
    The models are developed across multiple releases.
    The latest and most up to date reference is [TS38_901]_ section 7.3.
    The original development of the antenna models however, can be traced back to\
          [TS38_803]_ section 5.2.3.2.
    The antenna models are developed for deployment in different scenarios.

    Args:
        frequency (float): frequency of the antenna in [Hz]
        slant_angle (int): the polarization angle of the antenna [deg]
        polarization_model (int): polarization model accorsing to [TS38_901]_ section 7.3.2 [.]
        parameters_set (str): the parameter set describing the antenna power pattern. [.]
        resolution (:obj:`list` of two :obj:`int`): theta and phi steps for\
            pattern resolutions in degrees
        efficiency (float): total pattern efficiency. The pattern will be \
            normalized to this value. It is defined in percept from 0 to 100.
    """

    frequency: float
    slant_angle: int
    _slant_angle: float = field(init=False, repr=False)
    polarization_model: int
    _polarization_model: int = field(init=False, repr=False)
    parameters_set: str
    _parameters_set: str = field(init=False, repr=False)
    resolution: list[int]
    efficiency: float

    def __post_init__(self):
        self._pattern: PatternType = get_p_struct()

    # def _apply_antenna_polarization_model_1(
    #     self, power_pattern: npt.NDArray[np.float64]
    # ) -> tuple[npt.NDArray[np.complex64]]:
    #     sin_psi =

    #     e_th = np.sqrt(power_pattern) * np.cos(np.deg2rad(self.slant_angle))
    #     e_ph = np.sqrt(power_pattern) * np.sin(np.deg2rad(self.slant_angle))
    #     return e_th, e_ph

    def _apply_antenna_polarization_model_2(
        self, power_pattern: npt.NDArray[np.float64]
    ) -> tuple[npt.NDArray[np.complex64]]:
        e_th = np.sqrt(power_pattern) * np.cos(np.deg2rad(self.slant_angle))
        e_ph = np.sqrt(power_pattern) * np.sin(np.deg2rad(self.slant_angle))
        return e_th, e_ph

    def get_pattern(self) -> PatternType:
        """Generate analytical radiation pattern from 3GPP definitions.

        Generation is according to table 7.3-1 in [TS38_901]_.

        Returns:
            PatternType: pattern according to the parameters\
                :func:`~antenna_analysis.includes.utilities.get_p_struct`
        """

        thstep = self._get_th_step()
        phstep = self._get_th_step()
        vth, vph = self._get_vectors()

        vertical_cut = -np.minimum(
            12 * (((np.arange(0, 180 + thstep, thstep)) - 90) / self.theta_3db) ** 2,
            np.ones(len(vth)) * self.sla_v,
        )

        horizontal_cut = -np.minimum(
            (12 * (np.arange(-180, 180, phstep) / self.phi_3db) ** 2),
            np.ones(len(vph) - 1) * self.a_max,
        )

        v_patt = np.tile(vertical_cut, len(vph) - 1)
        v_patt = np.reshape(v_patt, (len(vth), len(vph) - 1), "F")
        h_patt = np.tile(horizontal_cut, len(vth))
        h_patt = np.reshape(h_patt, (len(vph) - 1, len(vth)), "F")
        h_patt = np.transpose(h_patt)

        power_patt = -np.minimum(
            -(v_patt + h_patt), np.ones((len(vth), len(vph) - 1)) * self.a_max
        )
        power_patt = 10 ** ((power_patt + self.g_e_max) / 10)

        if self.polarization_model == 1:
            e_th, e_ph = self._apply_antenna_polarization_model_2(power_patt)
            print(
                "Polarization model 1 requested but not implemented. Using polarization model 2."
            )
        elif self.polarization_model == 2:
            e_th, e_ph = self._apply_antenna_polarization_model_2(power_patt)
        else:
            sys.exit("unknown polarization model requested")

        # e_th, e_ph = self._normalize(e_th, e_ph)

        # Store the generated pattern
        self._pattern["source"] = "3GPP antenna element"
        self._pattern["file"] = "N/A"
        self._pattern["frequency"] = self.frequency
        self._pattern["port"] = 1
        self._pattern["th_step"] = self._resolution[0]
        self._pattern["ph_step"] = self._resolution[1]
        self._pattern["ff"] = np.stack((e_th, e_ph))

        return self._pattern

    @property
    def slant_angle(self) -> int:
        """Antenna polarization slant angle

        This is defined according to [TS38_901]_.
        The possible values are `0` and `45` degrees.
        `45` degrees correspond to a pair of cross-polarized antenna elements.
        `45` degrees really means `-45` for one antenna and `+45` degrees for the other.

        The default value is `0`.
        """
        return self._slant_angle

    @slant_angle.setter
    def slant_angle(self, value: int) -> None:
        self._slant_angle = value

    @property
    def polarization_model(self) -> int:
        """Set the polarization model of the antenna

        [TS38_901]_ defined two polarization models: 1 and 2.
        This parameter sets the model.

        The default value is model: 2.
        """
        return self._polarization_model

    @polarization_model.setter
    def polarization_model(self, value: int) -> None:
        self._polarization_model = value

    @property
    def parameters_set(self) -> str:
        """3GPP antenna parameter set.

        The 3GPP antenna patterns are defined by a set of parameters for different scenarios.
        Those are grouped here in tow main parameter sets: `outdoor` and `indoor`.

        `outdoor`(default) parameter set:
            
            The `outdoor` parameter set is the most commonly used antenna element pattern for \
                outdoor base stations. It can be traced back to [TS38_901]_ table 7.3-1,\
                [TS38_803]_ table 5.2.3.2.1-1 for urban macro scenario and table 5.2.3.2.2-1 for \
                dense urban scenario.

        `indoor` parameter set:

            The `indoor` parameter set is most commonly used antenna element pattern for indoor\
                base stations. It can be traced back to [TS38_803]_ table 5.2.3.2.3-1. \
                same parameters are also used in [TS38_803]_ table 5.2.3.3-1 as UE antenna pattern\
                definitions.
        """
        return self._parameters_set

    @parameters_set.setter
    def parameters_set(self, value: str) -> None:
        self._parameters_set = value

    @property
    def g_e_max(self) -> float:
        """Maximum antenna gain.

        For the `outdoor` parameter set this returns 8 dBi peak antenna gain.

        For the `indoor` parameter set this returns 5 dBi peak antenna gain.
        """
        if self.parameters_set == "outdoor":
            return 8  # dBi
        if self.parameters_set == "indoor":
            return 5  # dBi
        sys.exit("unknwon 3GPP parameter set")

    @property
    def theta_3db(self) -> float:
        """Vertical 3dB beamwidth.

        For the `outdoor` parameter set this returns 65 degrees beamwidth.

        For the `indoor` parameter set this returns 90 degrees beamwidth.
        """
        if self.parameters_set == "outdoor":
            return 65  # degrees
        if self.parameters_set == "indoor":
            return 90  # degrees
        sys.exit("unknwon 3GPP parameter set")

    @property
    def sla_v(self) -> float:
        """Vertical front to back ratio.

        For the `outdoor` parameter set this returns 30 dB.

        For the `indoor` parameter set this returns 25 dB.
        """
        if self.parameters_set == "outdoor":
            return 30  # dB
        if self.parameters_set == "indoor":
            return 25  # dB
        sys.exit("unknwon 3GPP parameter set")

    @property
    def phi_3db(self) -> float:
        """Horizontal 3dB beamwidth.

        For the `outdoor` parameter set this returns 65 degrees beamwidth.

        For the `indoor` parameter set this returns 90 degrees beamwidth.
        """
        if self.parameters_set == "outdoor":
            return 65  # degrees
        if self.parameters_set == "indoor":
            return 90  # degrees
        sys.exit("unknwon 3GPP parameter set")

    @property
    def a_max(self) -> float:
        """Horizontal front to back ratio.

        For the `outdoor` parameter set this returns 30 dB.

        For the `indoor` parameter set this returns 25 dB.
        """
        if self.parameters_set == "outdoor":
            return 30  # dB
        if self.parameters_set == "indoor":
            return 25  # dB
        sys.exit("unknwon 3GPP parameter set")


class TGPPAntennaBuilder:
    """Build a 3GPP antenna object."""

    def __init__(self):
        self._instance = None

    def __call__(
        self,
        frequency: float,
        slant_angle: int = 0,
        polarization_model: int = 2,
        parameters_set: str = "outdoor",
        resolution: list[int] = None,
        efficiency: float = 100,
        **_ignored,
    ) -> TGPPAntenna:
        """Build a 3GPP antenna object."""
        if resolution is None:
            resolution = [1, 1]
        self._instance = TGPPAntenna(
            frequency,
            slant_angle,
            polarization_model,
            parameters_set,
            resolution,
            efficiency,
        )
        return self._instance

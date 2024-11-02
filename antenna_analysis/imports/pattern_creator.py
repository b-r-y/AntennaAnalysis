"""Pattern creator module."""

from typing import Any
from antenna_analysis.imports.factory import ObjectFactory

from antenna_analysis.imports.pyramidal_horn import PyramidalHornBuilder
from antenna_analysis.imports.circular_patch import CircularPatchBuilder
from antenna_analysis.imports.tgpp_antenna import TGPPAntennaBuilder


class PatternCreator(ObjectFactory):
    """Pattern creator module."""

    def get(self, antenna_type: str, **kwargs: dict[str, Any]) -> Any:
        """Get the requested antenna type.

        Arguments:
            antenna_type (str): a string identifying the pattern to \
                be created. The possibilities are:

                * ``Pyramidal Horn`` - creates a pattern of type \
                    :class:`~antenna_analysis.imports.pyramidal_horn.PyramidalHorn`
                * ``Circular Patch`` - creates a pattern of type \
                    :class:`~antenna_analysis.imports.circular_patch.CircularPatch`
                * ``3GPP Antenna`` - creates a pattern of type \
                    :class:`~antenna_analysis.imports.tgpp_antenna.TGPPAntenna`


        Returns:
            pattern object  (PyramidalHorn, CircularPatch, TGPPAntenna): \
                an object of one of the possible patterns

        """
        return self.create(antenna_type, **kwargs)


pattern_creator = PatternCreator()

# Register antennas
pattern_creator.register_builder("Pyramidal Horn", PyramidalHornBuilder())
pattern_creator.register_builder("Circular Patch", CircularPatchBuilder())
pattern_creator.register_builder("3GPP Antenna", TGPPAntennaBuilder())

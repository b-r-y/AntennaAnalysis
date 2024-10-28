"""Pattern loaded from a file."""
from .common import CommonPattern
from ..includes.utilities import get_p_struct, PatternType


class FromFile(CommonPattern):
    """Define pattern loaded from a file."""

    def __init__(self, file: str, file_type: str):
        """Define pattern loaded from a file.

        Args:
            file (str): a string with the path to the source file
            file_type (str): identifier string for teh file type. Supported \
                types are:

                * `CST_File` - loads from a CST exported pattern. This type of\
                    input requires a file path.
                * `AAU_Satimo` - loads from AAU's Satimo chamber file exports.\
                    This type of input requires a file path.
        """
        self._file = file
        self._file_type = file_type
        self._pattern: PatternType = get_p_struct()

    def describe(self) -> str:
        """Return a description string.

        Generates and returns a description string of the antenna.

        Returns:
            str: string with summary of pattern properties.
        """
        string = (
            f'Pattern loaded from a file:\n'
            f'Source file: {self._file} \n'
            f'Source file type: {self._file_type} \n'
            f'Frequency: {self.frequency/1E6:.2f} [MHz]\n'
            f'Efficiency = {self.efficiency:.2f} [%]\n'
        )
        return string

    def get_pattern(self) -> PatternType:
        """Return the pattern structure.

        Returns:
            PatternType: pattern according to the structure\
                :func:`~antenna_analysis.includes.utilities.get_p_struct`
        """
        return self._pattern

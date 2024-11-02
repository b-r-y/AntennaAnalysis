from typing import Any, Union
from antenna_analysis.imports.pyramidal_horn import PyramidalHorn, PyramidalHornBuilder
from antenna_analysis.imports.circular_patch import CircularPatch, CircularPatchBuilder
from antenna_analysis.imports.tgpp_antenna import TGPPAntenna, TGPPAntennaBuilder


class ObjectFactory:
    def __init__(self):
        self._builders: dict[
            str, Union[PyramidalHornBuilder, CircularPatchBuilder, TGPPAntennaBuilder]
        ] = {}

    def register_builder(
        self,
        key: str,
        builder: Union[PyramidalHornBuilder, CircularPatchBuilder, TGPPAntennaBuilder],
    ):
        self._builders[key] = builder

    def create(
        self, key: str, **kwargs: Any
    ) -> Union[PyramidalHorn, CircularPatch, TGPPAntenna]:
        builder = self._builders.get(key)
        if not builder:
            raise ValueError(key)
        return builder(**kwargs)


factory = ObjectFactory()

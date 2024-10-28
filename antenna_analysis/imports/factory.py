from typing import Any, Union
from .pyramidal_horn import PyramidalHorn, PyramidalHornBuilder
from .circular_patch import CircularPatch, CircularPatchBuilder

class ObjectFactory:
    def __init__(self):
        self._builders: dict[str, Union[PyramidalHornBuilder,
                                        CircularPatchBuilder]] = {}

    def register_builder(self, key: str, builder: Union[PyramidalHornBuilder,
                                                        CircularPatchBuilder]):
        self._builders[key] = builder

    def create(self, key: str, **kwargs: Any) -> Union[PyramidalHorn,
                                                       CircularPatch]:
        builder = self._builders.get(key)
        if not builder:
            raise ValueError(key)
        return builder(**kwargs)

factory = ObjectFactory()

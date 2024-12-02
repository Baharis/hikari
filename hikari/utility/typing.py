import numbers
import os
import typing


PathLike = typing.Union[str, bytes, os.PathLike]
Number = typing.Union[typing.SupportsFloat, numbers.Number]

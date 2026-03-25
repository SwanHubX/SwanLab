from dataclasses import dataclass


@dataclass(frozen=True)
class SaveFile:
    source_path: str
    name: str
    target_path: str


__all__ = ["SaveFile"]

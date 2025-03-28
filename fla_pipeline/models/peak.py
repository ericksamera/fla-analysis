from dataclasses import dataclass
from typing import Optional


@dataclass
class Peak:
    position: float
    intensity: float
    saturated: bool = False
    note: Optional[str] = ""

    def to_dict(self) -> dict:
        return {
            "position": float(self.position),
            "intensity": float(self.intensity),
            "saturated": bool(self.saturated),
            "note": self.note or ""
        }

    @staticmethod
    def from_dict(data: dict) -> "Peak":
        return Peak(
            position=float(data["position"]),
            intensity=float(data["intensity"]),
            saturated=bool(data.get("saturated", False)),
            note=data.get("note", "")
        )

    @property
    def dict(self):
        return self.to_dict()

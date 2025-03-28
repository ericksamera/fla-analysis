# fla_pipeline/exporters/base.py

from abc import ABC, abstractmethod
import pandas as pd

class ExporterBase(ABC):
    @abstractmethod
    def export(self, df: pd.DataFrame) -> str:
        """Export to string format (CSV, TXT, etc.)"""
        pass

    @abstractmethod
    def filename(self) -> str:
        """Default filename for download"""
        pass

# fla_pipeline/core/fsa_loader.py

from typing import Dict
from Bio import SeqIO
import numpy as np

def load_fsa(file_path: str) -> Dict:
    """
    Loads an FSA file and extracts channel traces and size map.

    Args:
        file_path (str): Path to the .fsa or .ab1 file.

    Returns:
        Dict with keys:
            - name: str (filename)
            - smap: np.ndarray (size values)
            - channels: Dict[str, np.ndarray] (intensity values per dye channel)
    """
    fsa = SeqIO.read(file_path, 'abi')
    raw = fsa.annotations["abif_raw"]

    # Convert to numpy arrays
    smap = np.array(raw.get("SMap2", []), dtype=float)
    channels = {
        "6-FAM": np.array(raw.get("DATA9", []), dtype=float),
        "VIC":   np.array(raw.get("DATA10", []), dtype=float),
        "NED":   np.array(raw.get("DATA11", []), dtype=float),
        "PET":   np.array(raw.get("DATA12", []), dtype=float),
        "LIZ":   np.array(raw.get("DATA205", []), dtype=float),  # Size standard
    }

    return {
        "name": fsa.name,
        "smap": smap,
        "channels": channels
    }

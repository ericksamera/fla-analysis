# fla_pipeline/exporters/__init__.py

from .genalex import GenalexExporter

EXPORTERS = {
    "genalex": GenalexExporter,
    # "structure": StructureExporter,
    # "plink": PlinkExporter,
}

def get_exporter(name: str, **kwargs):
    exporter_cls = EXPORTERS.get(name.lower())
    if not exporter_cls:
        raise ValueError(f"Unknown exporter: {name}")
    return exporter_cls(**kwargs)
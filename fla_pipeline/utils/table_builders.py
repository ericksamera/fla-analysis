# fla_pipeline/utils/table_builders

import pandas as pd
from fla_pipeline.models.sample import Sample

def build_genotype_results_df(samples: dict[str, Sample]) -> pd.DataFrame:
    rows = []

    for sample in samples.values():
        sample_name = sample.metadata.get("Sample Name", sample.sample_id)
        for marker, genotype in sample.marker_results.items():
            rows.append({
                "sample_uid": sample.sample_uid,
                "sample_id": sample.sample_id,
                "Sample Name": sample_name,
                "Marker": marker,
                "Genotype": "/".join(map(str, genotype.alleles)),
                "Confidence": genotype.confidence,
                "QC Flags": "; ".join(genotype.qc_flags),
            })

    return pd.DataFrame(rows)

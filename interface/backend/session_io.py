# File: fla_pipeline/interface/backend/session_io.py

import json
import io
import streamlit as st
from typing import Dict
from fla_pipeline.models.sample import Sample
from fla_pipeline.models.genotype import GenotypeResult
from fla_pipeline.config.session_config import SessionConfig
from fla_pipeline.config import __VERSION__
from fla_pipeline.config.marker_config import MarkerConfig
from fla_pipeline.config.global_config import GlobalConfig
from fla_pipeline.analysis.config import DistanceConfig

def normalize_session_state():
    """Ensure key session_state entries are the correct object types."""

    if "marker_list" in st.session_state:
        st.session_state["marker_list"] = [
            MarkerConfig(**m) if isinstance(m, dict) else m
            for m in st.session_state["marker_list"]
        ]
    if "config" in st.session_state and isinstance(st.session_state["config"], dict):
        st.session_state["config"] = GlobalConfig(**st.session_state["config"])
    if "distance_config" in st.session_state and isinstance(st.session_state["distance_config"], dict):
        st.session_state["distance_config"] = DistanceConfig(**st.session_state["distance_config"])


def sanitize_for_json(obj):
    """Recursively sanitize any session object for JSON export."""
    if isinstance(obj, (int, float, str, bool)) or obj is None:
        return obj
    elif isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [sanitize_for_json(v) for v in obj]
    elif isinstance(obj, Sample):
        sample_dict = obj.to_dict()
        sample_dict["marker_results"] = {
            k: sanitize_for_json(v) for k, v in (obj.marker_results or {}).items()
        }
        return {"_type": "Sample", "data": sample_dict}
    elif isinstance(obj, GenotypeResult):
        return {"_type": "GenotypeResult", "data": obj.dict()}
    elif hasattr(obj, "model_dump"):
        return {"_type": obj.__class__.__name__, "data": obj.model_dump()}
    elif hasattr(obj, "dict"):
        return sanitize_for_json(obj.dict())
    return str(obj)


def deserialize_session(data: Dict):
    # Always ensure proper model types even if the loaded JSON has raw dicts
    session = SessionConfig(**data)

    st.session_state["samples"] = {
        s["sample_id"]: Sample.from_dict(s) for s in session.samples
    }

    st.session_state["marker_list"] = [
        m if isinstance(m, MarkerConfig) else MarkerConfig(**m)
        for m in session.marker_list
    ]

    st.session_state["config"] = (
        session.global_config
        if isinstance(session.global_config, GlobalConfig)
        else GlobalConfig(**session.global_config)
    )

    st.session_state["distance_config"] = session.distance_config

@st.cache_data(show_spinner="Serializing session…")
def serialize_session_cached(
    samples_json: list,
    marker_list_json: list,
    config_json: dict,
    distance_config_json: dict | None = None
) -> dict:
    session_config = SessionConfig(
        version=__VERSION__,
        samples=samples_json,
        marker_list=marker_list_json,
        global_config=config_json,
        distance_config=distance_config_json,
    )
    return session_config.model_dump()


def serialize_session() -> Dict:
    # Defensive guards to ensure session_state is clean
    if "config" not in st.session_state or not isinstance(st.session_state["config"], GlobalConfig):
        st.session_state["config"] = GlobalConfig()

    if "marker_list" not in st.session_state:
        st.session_state["marker_list"] = []

    if "samples" not in st.session_state:
        st.session_state["samples"] = {}

    # Now safe to pull values
    samples = st.session_state["samples"]
    marker_list = st.session_state["marker_list"]
    global_config = st.session_state["config"]
    distance_config = st.session_state.get("distance_config")

    session_config = SessionConfig(
        version=__VERSION__,
        samples=[s.to_dict() for s in samples.values()],
        marker_list=[
            m if isinstance(m, MarkerConfig) else MarkerConfig(**m)
            for m in marker_list
        ],
        global_config=global_config,
        distance_config=distance_config,
    )
    return session_config.model_dump()


def deserialize_session(data: Dict):
    """Load session from serialized session data dict."""
    session = SessionConfig(**data)
    st.session_state["samples"] = {
        s["sample_id"]: Sample.from_dict(s) for s in session.samples
    }
    st.session_state["marker_list"] = session.marker_list
    st.session_state["config"] = session.global_config
    st.session_state["distance_config"] = session.distance_config


# --- UI Hooks ---

def session_export_button():
    if st.button("Export Session", use_container_width=True, disabled=not bool(st.session_state.config)):

        samples_json = [s.to_dict() for s in st.session_state["samples"].values()]
        marker_list_json = [m.model_dump() for m in st.session_state["marker_list"]]
        config_json = st.session_state["config"].model_dump()
        distance_config_json = (
            st.session_state["distance_config"].model_dump()
            if st.session_state.get("distance_config") else None
        )

        session_dict = serialize_session_cached(
            samples_json,
            marker_list_json,
            config_json,
            distance_config_json
        )

        st.download_button(
            label="Download",
            data=json.dumps(session_dict, indent=2),
            file_name="fla_session.json",
            mime="application/json",
            use_container_width=True
        )



def session_import(file):
    """Read uploaded session JSON and populate session state."""
    try:
        raw = file.read()
        data = json.loads(raw.decode("utf-8"))
        deserialize_session(data)
        st.success("Session loaded.")
        st.rerun()
    except Exception as e:
        st.error(f"Failed to import session: {e}")


def session_import_button():
    """Show a button to open the session import dialog."""
    if st.button("Import Session", use_container_width=True):
        session_import_dialog()


@st.dialog("Import Session")
def session_import_dialog():
    """Dialog for uploading a session file."""
    uploaded = st.file_uploader("Upload session JSON", type="json")
    if uploaded:
        session_import(uploaded)

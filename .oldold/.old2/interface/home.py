#!/usr/bin/env python3
__description__ =\
"""
interface/home.py

Purpose:
  - Serves as the homepage for the FLA Viewer application.
  - Provides a brief introduction to the tool.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "Works"

import streamlit as st

def run():
    """Renders the home page."""
    st.title('Welcome to :rainbow: `abi-sauce`!')
    st.markdown('This tool processes pre-parsed FLA JSON files with extended raw data.')

run()
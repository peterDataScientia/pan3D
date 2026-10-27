# Pan3D — Streamlit Protein–Ligand 3D Figure Studio

Pan3D is a **Streamlit** web app for interactive protein–ligand binding site figures.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy with Streamlit Community Cloud

1. Visit https://share.streamlit.io/ and connect GitHub.
2. Select repository `peterDataScientia/pan3D`, branch `main`.
3. Set **Main file path** to `app.py` and deploy.

## Features in v0.2

- Upload complex PDB and optional spatially aligned ligand SDF/MOL/MOL2.
- Choose the HETATM ligand, visualize the 3D protein cartoon and ligand sticks.
- Detect residues within a selectable radius.
- Show pocket sticks and residue labels.
- Show **candidate N/O contacts** using distance only (NOT validated hydrogen bonds).
- Choose colors, rotate/zoom and export basic PNG from the embedded viewer.

## Limits

This is a first Streamlit prototype, not a fully validated publication-quality renderer. The dashed lines use proximity only; PLIP-backed contacts, draggable labels, ray-traced resolution-controlled images, PDF and SVG export are not implemented. A browser internet connection is required to load 3Dmol.js. Separate ligand files must have coordinates aligned with the protein.

The earlier static `index.html` remains as a legacy prototype; **`app.py` is the supported Streamlit entrypoint**.

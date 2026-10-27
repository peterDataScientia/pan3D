# Pan3D — Protein–Ligand 3D Figure Studio

Browser-based prototype for interactive binding-site illustrations.

## Current features
- Load a complex PDB or protein PDB with an aligned ligand in SDF, MOL2, or PDB format.
- Detect nonstandard HETATM ligands in a complex.
- Render protein cartoon, ligand sticks, nearby residue sticks and labels.
- Find residues within an adjustable distance cutoff.
- Show **candidate N/O proximity lines** at 2.2–3.5 Å (not verified hydrogen bonds).
- Configure colors, rotate/zoom and export a PNG.

## Run
Open `index.html` in a modern browser with an internet connection. The app loads 3Dmol.js from its CDN. Molecular files are read locally in the browser.

You may also serve the directory locally:
```sh
python -m http.server 8000
```
Then visit http://localhost:8000.

## Scientific limitations
- The first version is a **prototype**, not a validated molecular interaction profiler.
- Proximity lines do not evaluate donor/acceptor chemistry or angular geometry; do not label them as experimentally validated hydrogen bonds.
- When using a separate ligand file, its coordinates must be in the same frame as the protein.
- SVG/PDF export, publication-grade ray tracing, PLIP interaction classification, draggable residue labels and server-side PyMOL rendering are planned, not implemented.
- Test the rendering and exported graphics with real complexes before publication.

## Roadmap
1. PLIP-backed hydrogen bonds, hydrophobic contacts, salt bridges, pi interactions and validation.
2. Interactive label placement, interaction filtering and precise styling.
3. Publication-quality image export via PyMOL, plus vector/PDF outputs.
4. Quality assurance test set across complexes and input formats.

"""Pan3D: Streamlit protein–ligand interaction figure prototype."""
import json
from io import StringIO
import numpy as np
import streamlit as st
import streamlit.components.v1 as components
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

st.set_page_config(page_title="Pan3D | Protein–Ligand Figure Studio", page_icon="🧬", layout="wide")
st.title("Pan3D · Protein–Ligand Figure Studio")
st.caption("Streamlit molecular figure editor · 3Dmol.js renderer · v0.3")

AMINO = set("ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL MSE".split())
WATER = {"HOH", "WAT", "DOD"}
ELEMENTS = {"H", "D"}

def parse_pdb(data):
    structure = PDBParser(QUIET=True).get_structure("complex", StringIO(data))
    protein, ligands = [], {}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname.strip().upper()
                rid = f"{chain.id}:{residue.id[1]}{residue.id[2].strip()}:{resname}"
                atoms = []
                for atom in residue.get_atoms():
                    element = (atom.element or "").upper()
                    if element in ELEMENTS:
                        continue
                    atoms.append({"x": float(atom.coord[0]), "y": float(atom.coord[1]),
                                  "z": float(atom.coord[2]), "element": element,
                                  "name": atom.name, "residue": rid})
                if resname in AMINO:
                    protein.extend(atoms)
                elif resname not in WATER and residue.id[0].strip() and atoms:
                    ligands[rid] = atoms
        break
    return protein, ligands

def ligand_atoms_from_text(text, filename):
    # Parse with RDKit for SDF/MOL/MOL2 (coordinates must match protein frame).
    from rdkit import Chem
    lower = filename.lower()
    if lower.endswith(".mol2"):
        mol = Chem.MolFromMol2Block(text, removeHs=False, sanitize=False)
    else:
        mol = Chem.MolFromMolBlock(text, removeHs=False, sanitize=False)
    if mol is None or not mol.GetNumConformers():
        raise ValueError("Ligand file could not be parsed with 3D coordinates.")
    coords = mol.GetConformer().GetPositions()
    return [{"x": float(p[0]), "y": float(p[1]), "z": float(p[2]),
             "element": a.GetSymbol().upper(), "name": str(a.GetIdx()), "residue": "EXTERNAL"}
            for a, p in zip(mol.GetAtoms(), coords) if a.GetSymbol().upper() not in ELEMENTS]

def nearby(protein, ligand, cutoff):
    if not protein or not ligand:
        return [], []
    p = np.asarray([[a[k] for k in ("x", "y", "z")] for a in protein])
    l = np.asarray([[a[k] for k in ("x", "y", "z")] for a in ligand])
    tree = cKDTree(p)
    indices = set(i for point in l for i in tree.query_ball_point(point, cutoff))
    residues = sorted({protein[i]["residue"] for i in indices})
    # Proximity candidates, NOT chemically validated hydrogen bonds.
    contacts = []
    for j, a in enumerate(protein):
        if a["element"] not in ("N", "O"):
            continue
        for b in ligand:
            if b["element"] not in ("N", "O"):
                continue
            d = float(np.linalg.norm(np.array([a[k] for k in ("x", "y", "z")]) -
                                     np.array([b[k] for k in ("x", "y", "z")])))
            if 2.2 <= d <= 3.5:
                contacts.append({"start": {k: a[k] for k in ("x", "y", "z")},
                                 "end": {k: b[k] for k in ("x", "y", "z")}, "distance": round(d, 2)})
    return residues, contacts

with st.sidebar:
    st.header("1 · Molecular structures")
    pdb_file = st.file_uploader("Protein–ligand complex (PDB)", type=["pdb", "ent"])
    ligand_file = st.file_uploader("Optional aligned ligand (SDF/MOL/MOL2)", type=["sdf", "mol", "mol2"])
    st.header("2 · Binding pocket")
    cutoff = st.slider("Residue cutoff (Å)", 2.0, 10.0, 4.0, 0.5)
    show_residues = st.checkbox("Pocket residues as sticks", True)
    show_labels = st.checkbox("Residue labels", True)
    show_contacts = st.checkbox("Candidate N/O contacts", True)
    st.header("3 · Figure style")
    protein_color = st.color_picker("Protein cartoon", "#75BE87")
    ligand_color = st.color_picker("Ligand sticks", "#E66EBC")
    contact_color = st.color_picker("Contact dashes", "#D6B400")
    bg_color = st.color_picker("Background", "#FFFFFF")
    st.caption("Candidate contacts use only a distance check, not donor/acceptor geometry. Do not present them as validated H-bonds.")

if pdb_file is None:
    st.info("Upload a PDB complex (and optionally an aligned ligand) to generate a 3D binding-site illustration.")
    st.stop()

try:
    pdb_text = pdb_file.getvalue().decode("utf-8", errors="replace")
    protein, ligands = parse_pdb(pdb_text)
    extra_text, extra_ext = None, None
    if ligand_file is not None:
        extra_text = ligand_file.getvalue().decode("utf-8", errors="replace")
        extra_ext = ligand_file.name.rsplit(".", 1)[-1].lower()
        external_atoms = ligand_atoms_from_text(extra_text, ligand_file.name)
        ligands = {"External ligand": external_atoms, **ligands}
    if not ligands:
        st.error("No ligand HETATM residue detected. Upload an aligned SDF/MOL/MOL2 ligand.")
        st.stop()
except Exception as exc:
    st.error(f"Structure parsing failed: {exc}")
    st.stop()

ligand_name = st.selectbox("Select ligand", list(ligands))
ligand_atoms = ligands[ligand_name]
residues, contacts = nearby(protein, ligand_atoms, cutoff)
st.caption(f"{len(protein):,} protein heavy atoms · {len(ligand_atoms)} ligand heavy atoms · {len(residues)} pocket residues · {len(contacts)} candidate proximity lines")

# HTML rendering in an isolated Streamlit component. Molecular file contents are passed as JSON,
# not interpolated into JavaScript syntax.
payload = {
    "pdb": pdb_text,
    "ligandText": extra_text if ligand_name == "External ligand" else None,
    "ligandFormat": extra_ext,
    "selected": ligand_name,
    "residues": residues,
    "contacts": contacts if show_contacts else [],
    "showResidues": show_residues,
    "showLabels": show_labels,
    "proteinColor": protein_color,
    "ligandColor": ligand_color,
    "contactColor": contact_color,
    "background": bg_color,
    "ligandCoordinates": ligand_atoms,
}
blob = json.dumps(payload).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
html = """
<div style="font:12px sans-serif;color:#64748b;padding:4px">Drag to rotate · Scroll to zoom · Right drag to pan</div>
<div id="viewer" style="height:640px;width:100%;position:relative;background:#eef3f6;border:1px solid #cbd5e1"></div>
<div style="margin:9px 0"><button id="reset">Recenter</button> <button id="export">Export PNG</button> <span id="msg" style="margin-left:10px;color:#64748b"></span></div>
<script src="https://3Dmol.org/build/3Dmol-min.js" onerror="document.getElementById('msg').textContent='3Dmol CDN loading failed; check network';"></script>
<script>
const data = __DATA__;
const msg = document.getElementById('msg');
window.addEventListener('error', e => {msg.textContent='Pan3D rendering error: '+e.message;msg.style.color='crimson';});
if (typeof $3Dmol === 'undefined') { msg.textContent = '3Dmol.js failed to load. Check network or browser extensions.'; }
else {
  try {
  msg.textContent = 'Loading PDB atoms into 3D viewer...';
  const viewer = $3Dmol.createViewer(document.getElementById('viewer'), {backgroundColor: data.background, antialias:true});
  const model = viewer.addModel(data.pdb, 'pdb');
  if (!model.selectedAtoms({}).length) throw new Error('3Dmol could not read PDB atoms');
  model.setStyle({}, {cartoon:{color:data.proteinColor}});
  let ligandModel = model, selection = {};
  if (data.ligandText) {
    ligandModel = viewer.addModel(data.ligandText, data.ligandFormat === 'mol2' ? 'mol2' : 'sdf');
  } else {
    const parts = data.selected.split(':');
    // PDB style targets hetero residue by chain and residue number.
    selection = {chain:parts[0], resi:parseInt(parts[1], 10), resn:parts[2]};
  }
  if (data.showResidues) {
    for (const entry of data.residues) {
      const parts = entry.split(':');
      model.addStyle({chain:parts[0],resi:parseInt(parts[1],10),resn:parts[2]}, {stick:{radius:0.14}});
    }
  }
  if (data.ligandText) ligandModel.setStyle({}, {stick:{radius:0.22,color:data.ligandColor}});
  else model.addStyle(selection, {stick:{radius:0.22,color:data.ligandColor}});
  if (data.showLabels) {
    const atoms = model.selectedAtoms({});
    for (const entry of data.residues) {
      const parts = entry.split(':');
      const a = atoms.find(x=>x.chain===parts[0] && x.resi===parseInt(parts[1],10) && x.resn===parts[2] && x.atom==='CA');
      if(a) viewer.addLabel(parts[2]+parts[1], {position:{x:a.x,y:a.y,z:a.z},fontSize:13,fontColor:'#14213d',backgroundOpacity:0,inFront:true});
    }
  }
  for (const c of data.contacts) viewer.addCylinder({start:c.start,end:c.end,color:data.contactColor,radius:0.055,dashed:true});
  function center() {
    if (data.ligandText) viewer.zoomTo({model:ligandModel});
    else viewer.zoomTo(selection);
    viewer.zoom(0.75); viewer.render();
  }
  center();
  msg.textContent = 'Loaded: '+model.selectedAtoms({}).length+' atoms. Rotate the structure using your mouse.';
  document.getElementById('reset').onclick=center;
  document.getElementById('export').onclick=()=>{
    try {
      const a=document.createElement('a');
      a.href=viewer.pngURI();a.download='Pan3D_binding_site.png';
      document.body.appendChild(a);a.click();a.remove();
    }catch(e){msg.textContent='PNG export failed: '+e.message;}
  };
  } catch(error) { msg.textContent='Pan3D rendering error: '+error.message;msg.style.color='crimson';console.error(error); }
}
</script>
""".replace("__DATA__", blob)
components.html(html, height=730, scrolling=False)
with st.expander("Pocket residue identifiers"):
    st.write(", ".join(residues) if residues else "No residues found within the cutoff.")
with st.expander("Scientific interpretation"):
    st.warning("Yellow dashed lines are geometric N/O proximity candidates only. They are not PLIP-validated hydrogen bonds. Use validated interaction analysis before publication.")

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
st.caption("Streamlit molecular figure editor · 3Dmol.js renderer · v0.7")

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
    figure_mode = st.toggle("Publication figure mode", value=True)
    label_size = st.slider("Residue label size", 10, 24, 16)
    zoom_factor = st.slider("Camera zoom", 0.65, 1.4, 0.85, 0.05)
    cartoon_opacity = st.slider("Protein cartoon opacity", 0.1, 1.0, 0.82, 0.05)
    ligand_radius = st.slider("Ligand stick radius (Å)", 0.08, 0.30, 0.16, 0.01)
    pocket_radius = st.slider("Pocket residue stick radius (Å)", 0.06, 0.24, 0.10, 0.01)
    label_limit = st.slider("Maximum residue labels", 3, 35, 9)
    pocket_only = st.checkbox("Continuous protein cartoon (reference style)", True)
    rotation_x = st.slider("Camera tilt X (degrees)", -180, 180, 0, 10)
    rotation_y = st.slider("Camera turn Y (degrees)", -180, 180, 0, 10)
    rotation_z = st.slider("Camera roll Z (degrees)", -180, 180, 0, 10)
    preserve_elements = st.checkbox("Element colors for ligand O/N/S", True)
    label_offset = st.slider("Label outward offset (Å)", 0.0, 4.0, 1.5, 0.25)
    contact_radius = st.slider("Interaction dash thickness", 0.035, 0.14, 0.075, 0.005)
    export_scale = st.select_slider("PNG export scale", options=[1, 2, 3], value=2)
    st.caption("Final PNG is browser-rendered, not PyMOL ray-traced.")
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
# Rank labels by polar-contact proximity, then nearest atom distance.
lig_xyz = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms])
ranking = {}
for atom in protein:
    if atom["residue"] not in residues:
        continue
    pos = np.array([atom["x"], atom["y"], atom["z"]])
    distance = float(np.min(np.linalg.norm(lig_xyz - pos, axis=1)))
    ranking[atom["residue"]] = min(distance, ranking.get(atom["residue"], float("inf")))
polar_priority = set()
for atom in protein:
    if atom["residue"] in residues and atom["element"] in ("N", "O"):
        xyz = np.array([atom["x"], atom["y"], atom["z"]])
        for lig in ligand_atoms:
            if lig["element"] in ("N", "O"):
                dist = float(np.linalg.norm(xyz - [lig["x"], lig["y"], lig["z"]]))
                if 2.2 <= dist <= 3.5:
                    polar_priority.add(atom["residue"])
residues.sort(key=lambda r: (r not in polar_priority, ranking.get(r, 999), r))
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
    "figureMode": figure_mode,
    "labelSize": label_size,
    "zoomFactor": zoom_factor,
    "exportScale": export_scale,
    "cartoonOpacity": cartoon_opacity,
    "ligandRadius": ligand_radius,
    "pocketRadius": pocket_radius,
    "labelLimit": label_limit,
    "contactRadius": contact_radius,
    "pocketOnly": pocket_only,
    "rotationX": rotation_x,
    "rotationY": rotation_y,
    "rotationZ": rotation_z,
    "preserveElements": preserve_elements,
    "labelOffset": label_offset,
}
blob = json.dumps(payload).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
html = """
<style>
body{font-family:Arial,Helvetica,sans-serif;background:#f8fafc;color:#344456;margin:0}
#scene{width:100%;height:620px;position:relative;background:white;border:1px solid #d7e2e9;border-radius:8px}
#status{padding:9px;font-size:12px;color:#526579;white-space:pre-line}
button{margin:8px 8px 0 0;padding:10px 15px;border:0;border-radius:7px;background:#17649b;color:white;cursor:pointer}
</style>
<div id="scene"></div>
<button id="reset">Focus binding pocket</button><button id="export">Export PNG</button>
<div id="status">Loading figure engine…</div>
<script>
const D=__DATA__;
const el=document.getElementById('scene'),status=document.getElementById('status');
let viewer;
function say(t,bad=false){status.textContent=t;status.style.color=bad?'#b42332':'#526579';}
window.addEventListener('error', e=>say('Pan3D rendering error: '+e.message,true));
function run(){
try{
if(!window.$3Dmol)throw Error('3Dmol.js failed to load (CDN/network).');
viewer=$3Dmol.createViewer(el,{backgroundColor:D.background,antialias:true});
const protein=viewer.addModel(D.pdb,'pdb');
const atoms=protein.selectedAtoms({});
if(!atoms.length)throw Error('PDB could not be parsed by 3Dmol.js.');
// Keep a single continuous backbone ribbon: per-residue cartoon restyling caused detached green blades.
protein.setStyle({},{cartoon:{color:D.proteinColor,opacity:D.cartoonOpacity,style:'oval',arrows:true}});
const elementColors={O:'#ec3030',N:'#365ae0',S:'#d8b036',P:'#e6a346',F:'#61a944'};
const bits=String(D.selected).split(':');
const sel={chain:bits[0],resi:parseInt(bits[1],10),resn:bits[2]};
let ligand=protein;
if(D.ligandText){
  ligand=viewer.addModel(D.ligandText,D.ligandFormat==='mol2'?'mol2':'sdf');
  if(!ligand.selectedAtoms({}).length)throw Error('External ligand has no 3D atoms.');
  ligand.setStyle({},{stick:{radius:D.ligandRadius,color:D.ligandColor}});
  if(D.preserveElements)for(const [elem,color] of Object.entries(elementColors))ligand.setStyle({elem},{stick:{radius:D.ligandRadius,color}});
}else{
  if(!protein.selectedAtoms(sel).length)throw Error('Selected ligand not found: '+D.selected);
  // SetStyle on model is supported; addStyle belongs to the VIEWER.
  viewer.addStyle(sel,{stick:{radius:D.ligandRadius,color:D.ligandColor}});
  if(D.preserveElements)for(const [elem,color] of Object.entries(elementColors))viewer.addStyle({...sel,elem},{stick:{radius:D.ligandRadius,color}});
}
const labels=[];
for(const entry of D.residues.slice(0,D.labelLimit)){
  const [chain,num,resn]=entry.split(':');
  const resi=parseInt(num,10);
  const rsel={chain,resi,resn};
  if(D.showResidues)viewer.addStyle(rsel,{stick:{radius:D.pocketRadius}});
  if(D.showLabels){
    const a=atoms.find(t=>t.chain===chain && t.resi===resi && t.resn===resn && t.atom==='CA') ||
            atoms.find(t=>t.chain===chain && t.resi===resi && t.resn===resn);
    if(a){
      const cen=D.ligandCoordinates.reduce((s,p)=>({x:s.x+p.x,y:s.y+p.y,z:s.z+p.z}),{x:0,y:0,z:0});
      for(const axis of ['x','y','z'])cen[axis]/=Math.max(1,D.ligandCoordinates.length);
      const dx=a.x-cen.x,dy=a.y-cen.y,dz=a.z-cen.z,norm=Math.hypot(dx,dy,dz)||1;
      const at={x:a.x+D.labelOffset*dx/norm,y:a.y+D.labelOffset*dy/norm,z:a.z+D.labelOffset*dz/norm};
      labels.push({name:resn+num,at});
    }
  }
}
for(const label of labels)viewer.addLabel(label.name,{
position:{x:label.at.x,y:label.at.y,z:label.at.z},
font:'Arial',fontSize:D.labelSize,fontColor:'#101d28',bold:true,
backgroundOpacity:0,borderOpacity:0,inFront:true
});
for(const c of D.contacts){
  // Spheres + short cylinders make dashed proximity contacts visible in 3D and in PNG exports.
  const a=c.start,b=c.end,dx=b.x-a.x,dy=b.y-a.y,dz=b.z-a.z;
  const length=Math.sqrt(dx*dx+dy*dy+dz*dz);if(length<0.001)continue;
  const segments=Math.floor(length/0.35);
  for(let j=0;j<segments;j+=2){
    const start=j/segments,end=Math.min((j+1)/segments,1);
    const at=t=>({x:a.x+t*dx,y:a.y+t*dy,z:a.z+t*dz});
    viewer.addCylinder({start:at(start),end:at(end),radius:D.contactRadius,color:D.contactColor,fromCap:1,toCap:1});
  }
}
function focus(){
  // Fit ligand and its surrounding pocket using coordinate bounds instead of excessive ligand-only zoom.
  const ligandPoints=D.ligandCoordinates||[];
  const nearAtoms=atoms.filter(a=>D.residues.some(r=>{const [ch,num,name]=r.split(':');return a.chain===ch&&a.resi===parseInt(num,10)&&a.resn===name;}));
  const points=ligandPoints.concat(nearAtoms);
  if(points.length){
    const min=[Infinity,Infinity,Infinity],max=[-Infinity,-Infinity,-Infinity];
    points.forEach(p=>{[p.x,p.y,p.z].forEach((v,i)=>{min[i]=Math.min(min[i],v);max[i]=Math.max(max[i],v)})});
    const center={x:(min[0]+max[0])/2,y:(min[1]+max[1])/2,z:(min[2]+max[2])/2};
    const diameter=Math.max(...max.map((x,i)=>x-min[i]),12);
    // Standard 3Dmol zoomTo uses an atom selection and avoids invalid coordinate selectors.
    viewer.zoomTo(D.ligandText?{model:ligand.getID()}:sel);
    viewer.zoom(Math.min(1.4, D.zoomFactor));
  }else viewer.zoomTo();
  // Explicit camera orientation controls allow reproducible composition.
  if(D.rotationX)viewer.rotate(D.rotationX,'x');
  if(D.rotationY)viewer.rotate(D.rotationY,'y');
  if(D.rotationZ)viewer.rotate(D.rotationZ,'z');
  viewer.render();
}
focus();
document.getElementById('reset').onclick=focus;
document.getElementById('export').onclick=()=>{
try{
 // Temporarily upscale render target; restore viewer dimensions after export.
 const scale=D.exportScale||1;
 const oldW=el.clientWidth, oldH=el.clientHeight;
 if(scale>1){el.style.width=(oldW*scale)+'px';el.style.height=(oldH*scale)+'px';viewer.resize();viewer.render();}
 const url=viewer.pngURI();
 if(scale>1){el.style.width='100%';el.style.height='620px';viewer.resize();viewer.render();}
 const link=document.createElement('a');link.href=url;link.download='Pan3D_figure_'+scale+'x.png';
 document.body.appendChild(link);link.click();link.remove();
 say('PNG export generated ('+scale+'× canvas). Verify final figure resolution.');
}catch(e){el.style.width='100%';el.style.height='620px';viewer.resize();viewer.render();say('Export error: '+e.message,true);}
};
say('Rendered '+atoms.length+' atoms · '+D.residues.length+' pocket residues · '+D.contacts.length+' candidate polar contacts. Note: contacts are distance candidates, not PLIP validated.');
}catch(err){say('Pan3D rendering error: '+err.message,true);console.error(err);}
}
const script=document.createElement('script');
script.src='https://3Dmol.org/build/3Dmol-min.js';script.onload=run;
script.onerror=()=>{script.src='https://cdn.jsdelivr.net/npm/3dmol@2.5.3/build/3Dmol-min.js';script.onerror=()=>say('Unable to load rendering library from either CDN.',true);script.onload=run;};
document.head.appendChild(script);
</script>
""".replace("__DATA__", blob)
components.html(html, height=730, scrolling=False)
with st.expander("Pocket residue identifiers"):
    st.write(", ".join(residues) if residues else "No residues found within the cutoff.")
with st.expander("Scientific interpretation"):
    st.warning("Yellow dashed lines are geometric N/O proximity candidates only. They are not PLIP-validated hydrogen bonds. Use validated interaction analysis before publication.")

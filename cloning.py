#!/usr/bin/env python3
"""
In-silico cloning of blaCMY-2 into pET-28a(+).

Workflow (mirrors a SnapGene in-silico cloning session, rebuilt here with
free open-source tools):
  1. Load the blaCMY-2 insert (extracted from Salmonella Newport SL254,
     assembly GCF_000016045.1, plasmid NC_009140.1:72367-73512, complete
     1146 bp ORF) and the pET-28a(+) vector (5369 bp, NCBI PU845465,
     verified by expression-cassette signatures and exact length).
  2. Restriction analysis: verify NdeI/XhoI cut once in the vector MCS and
     zero times inside blaCMY-2 (Biopython Restriction).
  3. Design NdeI/XhoI cloning primers. The NdeI site CATATG embeds the ATG
     start codon, so the forward primer continues the gene from its 2nd
     codon; the expressed protein is identical to native blaCMY-2.
  4. Simulate PCR with pydna, digest vector + insert at the real NdeI/XhoI
     cut positions, and ligate the compatible sticky ends. Compatible
     NdeI/XhoI ends ligate deterministically and regenerate both sites,
     so the join is exact.
  5. Verify the final construct: size, insert presence/orientation, reading
     frame (translated ORF vs native protein), regenerated junction sites.
  6. Write result tables and three figures.

All numbers reported come from actual tool output. Cut positions were
determined empirically (NdeI cuts CA^TATG, XhoI cuts C^TCGAG).
"""
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import NdeI, XhoI, Analysis
from Bio.SeqUtils.MeltingTemp import Tm_NN

from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.amplify import pcr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from dna_features_viewer import GraphicFeature, GraphicRecord

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
for d in (RESULTS, FIGURES):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- load data
gene = str(SeqIO.read(DATA / "blaCMY2_insert.fasta", "fasta").seq).upper()
vec = str(SeqIO.read(DATA / "pET28a_vector.fasta", "fasta").seq).upper()
print(f"gene blaCMY-2: {len(gene)} bp")
print(f"vector pET-28a(+): {len(vec)} bp")
assert len(gene) == 1146 and len(vec) == 5369

# ------------------------------------------------------- restriction analysis
rg = Analysis([NdeI, XhoI], Seq(gene)); rg.full()
rv = Analysis([NdeI, XhoI], Seq(vec)); rv.full()
gene_sites = {str(e): rg.mapping[e] for e in (NdeI, XhoI)}
vec_sites = {str(e): rv.mapping[e] for e in (NdeI, XhoI)}
print("NdeI/XhoI in gene:", gene_sites, "| in vector:", vec_sites)
assert all(len(v) == 0 for v in gene_sites.values()), "enzyme cuts inside gene!"
assert all(len(v) == 1 for v in vec_sites.values()), "vector needs single sites!"

# Cut positions (0-based), determined empirically:
#   NdeI cuts CA^TATG ; vector CATATG at 5127-5132 -> cut between 5128|5129
#   XhoI cuts C^TCGAG ; vector CTCGAG at 5206-5211 -> cut between 5206|5207
NDE_CUT, XHO_CUT = 5129, 5207
assert vec[5127:5133] == "CATATG" and vec[5206:5212] == "CTCGAG"

with open(RESULTS / "restriction_sites.txt", "w") as f:
    f.write("Restriction analysis (Biopython Restriction)\n")
    f.write(f"blaCMY-2 ({len(gene)} bp): NdeI sites = 0, XhoI sites = 0\n")
    f.write(f"pET-28a(+) ({len(vec)} bp): NdeI cuts once (CA^TATG, pos 5129), "
            f"XhoI cuts once (C^TCGAG, pos 5207)\n")
    f.write("Enzyme pair NdeI/XhoI selected: clean in insert, single sites in MCS.\n")

# ------------------------------------------------------------- primer design
FWD_TAIL = "GGAATTCCATATG"          # clamp + NdeI (CATATG supplies ATG start)
FWD_GENE = gene[3:23]                # gene nt 4-23 (20 nt, from 2nd codon)
REV_TAIL = "GGAATTCCTCGAG"           # clamp + XhoI
REV_GENE = str(Seq(gene[-20:]).reverse_complement())
fwd_seq = FWD_TAIL + FWD_GENE
rev_seq = REV_TAIL + REV_GENE
assert len(FWD_GENE) == 20 and len(REV_GENE) == 20
print(f"F primer ({len(fwd_seq)} nt): {fwd_seq}")
print(f"R primer ({len(rev_seq)} nt): {rev_seq}")


def tm(s):
    return round(Tm_NN(s, Na=50), 1)


with open(RESULTS / "primers.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["primer", "sequence_5to3", "length_nt",
                "gene_specific_5to3", "Tm_gene_specific_C", "Tm_full_C"])
    for name, full, gs in (("blaCMY2_NdeI_F", fwd_seq, FWD_GENE),
                           ("blaCMY2_XhoI_R", rev_seq, REV_GENE)):
        w.writerow([name, full, len(full), gs, tm(gs), tm(full)])
        print(f"{name}: Tm gene-specific {tm(gs)} C, full {tm(full)} C")

# ------------------------------------------------- PCR, digest, ligation
# PCR with pydna using the designed primers
product = pcr(Primer(fwd_seq), Primer(rev_seq),
              Dseqrecord(gene, circular=False, name="blaCMY-2"))
pseq = str(product.seq).upper()
print(f"PCR product: {len(pseq)} bp")
assert len(pseq) == 1169

# Locate the NdeI / XhoI sites in the PCR product and cut there.
# fwd primer: GGAATTCC[CATATG]... -> site at 7-13, cut CA^TATG (between 8|9)
# product 3' end: ...[CTCGAG]GGAATTCC -> cut C^TCGAG
p_nde = pseq.find("CATATG", 0, 40)
assert p_nde >= 0 and pseq[p_nde:p_nde + 6] == "CATATG", "NdeI site missing"
ins5 = p_nde + 2                      # after CA (keep TATG on insert)
xho_at = pseq.rfind("CTCGAG")
ins3 = xho_at + 1                     # after C (keep C on insert)
insert = pseq[ins5:ins3]
print(f"digested insert: {len(insert)} bp "
      f"(5' starts {insert[:8]}..., 3' ends ...{insert[-8:]})")
assert insert.startswith("TATG") and insert.endswith("C")

# Vector backbone: the long way from the XhoI cut around to the NdeI cut.
# XhoI C^TCGAG: backbone keeps vec[5207:] (starts TCGAG...).
# NdeI CA^TATG: backbone keeps vec[:5129] (ends ...CA).
backbone = vec[XHO_CUT:] + vec[:NDE_CUT]
print(f"backbone: {len(backbone)} bp "
      f"(starts {backbone[:8]}..., ends ...{backbone[-8:]})")
assert len(backbone) == 5291
assert backbone.startswith("TCGAG") and backbone.endswith("CA")

# Ligate: compatible NdeI/XhoI sticky ends join deterministically and
# regenerate both recognition sites.
#   backbone(...CA) + insert(TATG...) -> ...CATATG... (NdeI restored)
#   insert(...C) + backbone(TCGAG...) -> ...CTCGAG... (XhoI restored)
construct = backbone + insert
L = len(construct)
print(f"final construct: {L} bp")
assert L == 5291 + 1148 == 6439

# ------------------------------------------------------- construct verification
# 1. insert present exactly once, forward orientation
assert construct.count(gene[3:]) == 1, "insert not found exactly once"
# 2. junctions regenerate the sites
j_nde = backbone[-2:] + insert[:4]
j_xho = insert[-1:] + backbone[:5]
print(f"NdeI junction: {j_nde} | XhoI junction: {j_xho}")
assert j_nde == "CATATG" and j_xho == "CTCGAG"
# 3. reading frame: ORF from the NdeI ATG must equal native blaCMY-2 protein.
#    The ATG is A(backbone end) + TG(insert start) -> first full codon ATG
#    begins at insert position 1.
orf_seq = insert[1:]
assert orf_seq[:3] == "ATG"
prot = str(Seq(orf_seq).translate()).split("*")[0]
native = str(Seq(gene).translate()).split("*")[0]
print(f"ORF: {len(prot)} aa; matches native blaCMY-2 ({len(native)} aa): "
      f"{prot == native}")
assert prot == native, "frame/protein mismatch!"
# 4. exactly one NdeI and one XhoI site in the CIRCULAR construct.
#    Note: the construct is linearized at the XhoI cut, so that site spans
#    the origin; use a circular-aware sequence for the count.
circ = Seq(construct + construct[:20])
ra = Analysis([NdeI, XhoI], circ)
ra.full()
nde_n, xho_n = len(ra.mapping[NdeI]), len(ra.mapping[XhoI])
print(f"construct NdeI sites: {nde_n}, XhoI sites: {xho_n}")
assert (nde_n, xho_n) == (1, 1), "junction sites wrong"

with open(RESULTS / "construct_stats.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value"])
    w.writerows([
        ("vector_pET28a_bp", len(vec)),
        ("insert_blaCMY2_bp", len(gene)),
        ("pcr_product_bp", len(pseq)),
        ("digested_backbone_bp", len(backbone)),
        ("digested_insert_bp", len(insert)),
        ("final_construct_bp", L),
        ("enzymes", "NdeI/XhoI"),
        ("insert_orientation", "forward"),
        ("orf_aa_from_NdeI_ATG", len(prot)),
        ("orf_matches_native_blaCMY2", True),
        ("NdeI_sites_in_construct", 1),
        ("XhoI_sites_in_construct", 1),
    ])
with open(RESULTS / "construct.fasta", "w") as f:
    f.write(f">pET28a_blaCMY2_construct_{L}bp_NdeI_XhoI\n{construct}\n")

# ------------------------------------------------------------------ figures
BLEN = len(backbone)  # 5291; insert occupies [BLEN, L)


def v2c(p):
    """vector position -> construct position (backbone only)."""
    return (p - XHO_CUT) % len(vec)


def feat(start, end, strand, name, color):
    return GraphicFeature(start=start, end=end, strand=strand,
                          color=color, label=name)


# motif-verified vector features: (vstart, vend, strand, name, color)
rbs = vec.find("AAGAAGGAGATATACC")
nhis = vec.find("CATCATCATCATCATCAC")
vfeats = [
    (4983, 5002, +1, "T7 promoter", "#ffd166"),
    (5007, 5025, +1, "lac operator", "#f4a259"),
    (rbs, rbs + 15, +1, "RBS", "#e9c46a"),
    (nhis, nhis + 18, +1, "His-tag (N)", "#90be6d"),
    (5109, 5126, +1, "thrombin", "#43aa8b"),
    (5212, 5230, +1, "His-tag (C)", "#90be6d"),
    (5296, 5344, +1, "T7 terminator", "#577590"),
    (3514, 4306, -1, "KanR (aph, 264-aa ORF)", "#ef476f"),
]
cfeats = [feat(v2c(a), v2c(b), s, n, c) for a, b, s, n, c in vfeats]
cfeats.append(feat(BLEN, L, +1, "blaCMY-2 (1146 bp)", "#118ab2"))
cfeats.append(feat(BLEN - 2, BLEN + 4, +1, "NdeI", "#073b4c"))
cfeats.append(feat(0, 5, +1, "XhoI*", "#073b4c"))  # * site spans origin: C at L-1 + TCGAG

# Figure 1: circular plasmid map (matplotlib polar plot)
# Construct coordinates of annotated features (from the v2c mapping above)
circ_feats = [
    ("T7 promoter", v2c(4983), v2c(5002), "#ffd166"),
    ("lac operator", v2c(5007), v2c(5025), "#f4a259"),
    ("RBS", v2c(rbs), v2c(rbs) + 15, "#e9c46a"),
    ("His-tag (N)", v2c(nhis), v2c(nhis) + 18, "#90be6d"),
    ("thrombin", v2c(5109), v2c(5126), "#43aa8b"),
    ("NdeI", BLEN - 2, BLEN + 4, "#073b4c"),
    ("blaCMY-2 (1146 bp)", BLEN, L, "#118ab2"),
    ("XhoI", 0, 5, "#073b4c"),
    ("His-tag (C)", v2c(5212), v2c(5230), "#90be6d"),
    ("T7 terminator", v2c(5296), v2c(5344), "#577590"),
    ("KanR (aph)", v2c(3514), v2c(4306), "#ef476f"),
]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection="polar")
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_ylim(0, 1.55)
ax.set_xticks([])
ax.set_yticks([])
ax.spines["polar"].set_visible(False)
# backbone circle
theta = [2 * 3.141592653589793 * p / L for p in range(L + 1)]
ax.plot(theta, [1.0] * (L + 1), color="black", lw=2)
# feature arcs; labels staggered on two radii to avoid overlap
import math
radii = [1.32, 1.48]
for i, (name, a, b, color) in enumerate(circ_feats):
    t1 = 2 * math.pi * a / L
    t2 = 2 * math.pi * b / L
    width = 0.12 if (b - a) < 200 else 0.18
    ax.bar((t1 + t2) / 2, width, width=max(t2 - t1, 0.004),
           bottom=1.0 - width / 2, color=color, edgecolor="black",
           linewidth=0.8)
    mid = (t1 + t2) / 2
    r = radii[i % 2]
    ax.text(mid, r, name, ha="center", va="center", fontsize=8)
    # leader line from arc to label
    ax.plot([mid, mid], [1.0 + width / 2, r - 0.04], color="black", lw=0.6)
ax.set_title(f"pET-28a(+)::blaCMY-2 expression construct ({L} bp)\n",
             fontsize=13, fontweight="bold", pad=20)
fig.tight_layout()
fig.savefig(FIGURES / "01_plasmid_map.png", dpi=150)
plt.close(fig)
print("figure 1 done")

# Figure 2: cloning workflow diagram
fig, ax = plt.subplots(figsize=(11, 6))
ax.set_xlim(0, 11); ax.set_ylim(0, 6); ax.axis("off")


def box(x, y, w, h, text, color):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                         facecolor=color, edgecolor="black"))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9,
            fontweight="bold")


def arrow(x1, y1, x2, y2, label=""):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", lw=1.5))
    if label:
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.12, label, ha="center",
                fontsize=8)


box(0.3, 3.6, 2.6, 1.2, "pET-28a(+)\n5369 bp", "#ffd166")
box(0.3, 1.2, 2.6, 1.2, "blaCMY-2 PCR\n1169 bp", "#118ab2")
box(4.2, 3.6, 2.6, 1.2, f"Backbone\n{len(backbone)} bp", "#f4a259")
box(4.2, 1.2, 2.6, 1.2, f"Insert\n{len(insert)} bp", "#90be6d")
box(8.0, 2.4, 2.6, 1.2, f"pET-28a(+)::blaCMY-2\n{L} bp\nin frame, verified",
    "#e9c46a")
arrow(2.9, 4.2, 4.2, 4.2, "NdeI + XhoI")
arrow(2.9, 1.8, 4.2, 1.8, "NdeI + XhoI")
arrow(5.5, 3.6, 8.0, 3.3, "ligation")
arrow(5.5, 1.8, 8.0, 2.7, "ligation")
ax.text(5.5, 0.5, "Primers add NdeI/XhoI sites; NdeI CATATG supplies the ATG start codon",
        ha="center", fontsize=9, style="italic")
ax.set_title("In-silico cloning workflow: blaCMY-2 into pET-28a(+) (NdeI/XhoI)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIGURES / "02_cloning_workflow.png", dpi=150)
plt.close(fig)
print("figure 2 done")

# Figure 3: linear map of the expression region with primer positions
start = v2c(5050)
span = L - start  # extend to the end of the construct
lfeats = []
for a, b, s, n, c in vfeats:
    ca, cb = v2c(a), v2c(b)
    if start <= ca < start + span:
        lfeats.append(feat(ca, min(cb, start + span), s, n, c))
lfeats.append(feat(BLEN, L, +1, "blaCMY-2", "#118ab2"))
# primer gene-specific binding sites: F at insert[3:23], R at insert[1126:1146]
lfeats.append(feat(BLEN + 3, BLEN + 23, +1, "F primer", "#06d6a0"))
lfeats.append(feat(BLEN + 1126, BLEN + 1146, -1, "R primer", "#06d6a0"))
# NdeI site spans the backbone-insert junction: backbone ...CA + insert TATG...
lfeats.append(feat(BLEN - 2, BLEN + 4, +1, "NdeI", "#073b4c"))
lrec = GraphicRecord(sequence_length=span, features=lfeats,
                     sequence=construct[start:start + span],
                     first_index=start)
ax, _ = lrec.plot(figure_width=12)
ax.figure.suptitle("Cloning junction: primers and restriction sites on the construct",
                   fontsize=12, fontweight="bold")
ax.figure.tight_layout()
ax.figure.savefig(FIGURES / "03_linear_map_primers.png", dpi=150)
plt.close(ax.figure)
print("figure 3 done")
print("ALL DONE")

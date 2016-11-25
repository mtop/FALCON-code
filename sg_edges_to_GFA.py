#!/usr/local/packages/anaconda2/bin/python
## This is the `sg_edges_to_GFA.py` script
## (More) information at https://github.com/PacificBiosciences/FALCON/wiki/Convert-FALCON-assembly-graph-to-GFA-format

from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader

read_in_graph = set()
G_asm = AsmGraph("sg_edges_list", "utg_data", "ctg_paths")
edge_to_ctg = {}
a_path = {}
with open("a_ctg_tiling_path") as f:
    for row in f:
        row = row.strip().split()
        ctg_id, v, w, edge_rid, b, e  = row[:6]
        a_path.setdefault( ctg_id, [] )
        a_path[ctg_id].append( (v, w) )
        ctg_id = ctg_id.split("-")[0] #get the primary contig id
        edge_to_ctg[ (v, w) ] = ctg_id, "A"

p_path = {}
with open("p_ctg_tiling_path") as f:
    for row in f:
        row = row.strip().split()
        ctg_id, v, w, edge_rid, b, e  = row[:6]
        p_path.setdefault( ctg_id, [] )
        p_path[ctg_id].append( (v, w) )
        edge_to_ctg[ (v, w) ] = ctg_id, "P"

read_pairs = set()
link_lines = []
for v, w in G_asm.sg_edges:
    edge_data = G_asm.sg_edges[(v, w)]
    """000084959:E 000376804:B 000376804 25867 0 1766 99.55 TR"""
    if edge_data[-1] == "G":
        r1, r1end = v.split(":")
        r2, r2end = w.split(":")
        rp = [r1, r2]
        rp.sort()
        rp = tuple(rp)
        if rp in read_pairs:  #GFA model read rather than the end, so we need to collapse dual edge
            continue
        else:
            read_pairs.add( rp )
        read_in_graph.add(r1)
        read_in_graph.add(r2)
        if r1end == "E":
            o1 = "+"
        else:
            o1 = "-"
        if r2end == "E":
            o2 = "+"
        else:
            o2 = "-"
        overlap_length = int(edge_data[1])
        overlap_idt = float(edge_data[2])
        ctg_id = edge_to_ctg.get( (v, w), ("NA","NA") )
        link_lines.append( "\t".join( ["L", r1, o1, r2, o2, "*",
            "ol:i:%d" % overlap_length, "oi:f:%.1f" % overlap_idt,
            "ci:A:%s-%s" % ctg_id]  ) )

try:
	# Works with v1.7.5
	f = FastaReader("../1-preads_ovl/db2falcon/preads4falcon.fasta")
except:
	try:
		# Works with v1.8.2
		f = FastaReader("../1-preads_ovl/preads4falcon.fasta")
	except:
		print "Bummer, this code does not work with your version of FALCON."

seq_len = {}
for r in f:
    if r.name not in read_in_graph:
        continue
    seq_len[r.name] = len(r.sequence)

print "H\tVN:Z:1.0"
read_in_graph = list(read_in_graph)
read_in_graph.sort()
for r in read_in_graph:
    print "\t".join( ["S", r, "*", "LN:i:%s" % seq_len[r]] )
for link_line in link_lines:
    print link_line

p_path_k = p_path.keys()
p_path_k.sort()
for k in p_path_k:
    v, w = p_path[k][0]
    rn, rend = v.split(":")
    o = "+" if rend == "E" else "-"
    segs = [rn+o]
    for v, w in p_path[k]:
        rn, rend = w.split(":")
        o = "+" if rend == "E" else "-"
        segs.append( rn+o )
    out = ["P", k, ",".join(segs), ",".join(["*"]*len(segs))] #yes, I know it is a little bit ....
    print "\t".join(out)

a_path_k = a_path.keys()
a_path_k.sort()
for k in a_path_k:
    v, w = a_path[k][0]
    rn, rend = v.split(":")
    o = "+" if rend == "E" else "-"
    segs = [rn+o]
    for v, w in a_path[k]:
        rn, rend = w.split(":")
        o = "+" if rend == "E" else "-"
        segs.append( rn+o )
    out = ["A", k, ",".join(segs), ",".join(["*"]*len(segs))] #yes, I know it is a little bit ....
    print "\t".join(out)

#!/usr/bin/env python3
"""Convertit un maillage Gmsh MSH 2.2 ASCII en maillage au format natif Aster (.mail)."""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ELEMENT_TYPES = {1: ("SEG2", 2), 2: ("TRIA3", 3), 8: ("SEG3", 3), 9: ("TRIA6", 6), 15: ("POI1", 1), 16: ("QUAD8", 8)}


def _section(lines: list[str], name: str) -> list[str]:
    start = lines.index(f"${name}") + 1
    return lines[start : lines.index(f"$End{name}", start)]


def read_msh22(text: str):
    lines = [line.strip() for line in text.splitlines()]
    version = _section(lines, "MeshFormat")[0].split()
    if not version[0].startswith("2") or version[1] != "0":
        raise ValueError(f"MSH 2.x ASCII attendu, recu : {version}")
    names = {}
    for row in _section(lines, "PhysicalNames")[1:]:
        _dim, tag, name = row.split(maxsplit=2)
        names[int(tag)] = name.strip('"')
    nodes = {}
    for row in _section(lines, "Nodes")[1:]:
        tag, x, y, _z = row.split()
        nodes[int(tag)] = (float(x), float(y))
    elements, groups = [], defaultdict(list)
    for row in _section(lines, "Elements")[1:]:
        fields = [int(v) for v in row.split()]
        tag, kind, ntags = fields[:3]
        if kind not in ELEMENT_TYPES:
            raise ValueError(f"Type d'element Gmsh non gere : {kind}")
        phys = fields[3]
        connectivity = fields[3 + ntags :]
        if len(connectivity) != ELEMENT_TYPES[kind][1]:
            raise ValueError(f"Connectivite invalide pour l'element {tag}")
        elements.append((tag, ELEMENT_TYPES[kind][0], connectivity))
        groups[names[phys]].append(tag)
    return nodes, elements, groups


def write_aster(nodes, elements, groups) -> str:
    out = ["TITRE", "Eprouvette CT - genere par msh_to_aster.py", "FINSF", "COOR_2D"]
    out += [f"N{tag} {x:.15e} {y:.15e}" for tag, (x, y) in sorted(nodes.items())]
    out.append("FINSF")
    by_type = defaultdict(list)
    point_nodes = {}
    for tag, kind, connectivity in elements:
        if kind == "POI1":  # un point physique devient un groupe de noeuds, sans maille
            point_nodes[tag] = connectivity[0]
        else:
            by_type[kind].append((tag, connectivity))
    for kind, members in by_type.items():
        out.append(kind)
        out += [f"M{tag} " + " ".join(f"N{n}" for n in conn) for tag, conn in members]
        out.append("FINSF")
    for name, tags in groups.items():
        if all(tag in point_nodes for tag in tags):
            out.append(f"GROUP_NO NOM={name}")
            out += [f"N{point_nodes[tag]}" for tag in tags]
        else:
            out.append(f"GROUP_MA NOM={name}")
            out += [f"M{tag}" for tag in tags]
        out.append("FINSF")
    out.append("FIN")
    return "\n".join(out) + "\n"


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("usage : msh_to_aster.py entree.msh sortie.mail")
    nodes, elements, groups = read_msh22(Path(sys.argv[1]).read_text())
    Path(sys.argv[2]).write_text(write_aster(nodes, elements, groups))


if __name__ == "__main__":
    main()

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from msh_to_aster import read_msh22, write_aster
from run_etape1 import parse_table

MSH = """$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
3
0 1 "P"
1 2 "L"
2 3 "S"
$EndPhysicalNames
$Nodes
6
1 0 0 0
2 1 0 0
3 0 1 0
4 0.5 0 0
5 0.5 0.5 0
6 0 0.5 0
$EndNodes
$Elements
3
1 15 2 1 1 1
2 8 2 2 1 1 2 4
3 9 2 3 1 1 2 3 4 5 6
$EndElements
"""


def test_conversion_groupes():
    text = write_aster(*read_msh22(MSH))
    assert "GROUP_NO NOM=P\nN1\nFINSF" in text
    assert "GROUP_MA NOM=L\nM2\nFINSF" in text
    assert "TRIA6\nM3 N1 N2 N3 N4 N5 N6\nFINSF" in text
    assert "POI1" not in text


def test_parse_table_plusieurs_tables():
    text = "#\n#---\n##TITRE\nA B\n1 2.0E+00\n\n#\n##TITRE\nA B\n3 4.0E+00\n"
    assert parse_table(text) == [{"A": "1", "B": "2.0E+00"}, {"A": "3", "B": "4.0E+00"}]

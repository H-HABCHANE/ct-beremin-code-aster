// Demi-eprouvette CT (y >= 0), deformations planes. Unites : mm.
// Geometrie proportionnelle a W (ASTM E399) ; origine sur la ligne de charge.
DefineConstant[ W = 50, aW = 0.5, h_tip = 0.05, h_far = 2.5, R_fine = 1.0, quad = 0 ];
a = aW*W;
R = 0.125*W;      // rayon du trou de goupille
yc = 0.275*W;     // ordonnee du centre du trou
xf = -0.25*W;     // face avant
H = 0.6*W;        // demi-hauteur

Point(1) = {xf, 0, 0};  Point(2) = {a, 0, 0};  Point(3) = {W, 0, 0};
Point(4) = {W, H, 0};   Point(5) = {xf, H, 0};
Point(6) = {0, yc, 0};  // centre du trou (non maille)
Point(7) = {R, yc, 0};  Point(8) = {0, yc+R, 0};
Point(9) = {-R, yc, 0}; Point(10) = {0, yc-R, 0};

Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 5}; Line(5) = {5, 1};
Circle(6) = {7, 6, 8}; Circle(7) = {8, 6, 9}; Circle(8) = {9, 6, 10}; Circle(9) = {10, 6, 7};
Curve Loop(1) = {1, 2, 3, 4, 5};
Curve Loop(2) = {6, 7, 8, 9};
Plane Surface(1) = {1, 2};

Physical Surface("SOLIDE") = {1};
Physical Curve("LEVRE") = {1};
Physical Curve("LIGAMENT") = {2};
Physical Curve("TROU_HAUT") = {6, 7};
Physical Point("POINTE") = {2};
Physical Point("BLOC_X") = {3};
Physical Point("SOMMET_TROU") = {8};

Field[1] = Distance; Field[1].PointsList = {2};
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = h_tip; Field[2].SizeMax = h_far;
Field[2].DistMin = R_fine; Field[2].DistMax = 0.4*W;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6;
Mesh.ElementOrder = 2;
Mesh.MshFileVersion = 2.2;
If (quad)
  Recombine Surface {1};
  Mesh.RecombinationAlgorithm = 3;
  Mesh.SecondOrderIncomplete = 1;
EndIf

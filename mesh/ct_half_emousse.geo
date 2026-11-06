// Demi-eprouvette CT a pointe emoussee (rayon r0), mailles QUAD8. Unites : mm.
// Autour de la pointe : maillage structure en deux blocs (A : couronne, B : au-dessus de la levre).
DefineConstant[ W = 50, aW = 0.5, r0 = 0.001, Rf = 1.0, Lb = 0.6,
                n_theta = 32, n_r = 56, n_b = 28, q_r = 1.12, q_b = 1.17, h_far = 2.5 ];
a = aW*W;
R = 0.125*W; yc = 0.275*W; xf = -0.25*W; H = 0.6*W;
xc = a - r0;  // centre de l'arrondi

// Pointe
Point(1) = {xc, 0, 0};            // centre de l'arrondi (non maille)
Point(2) = {a, 0, 0};             // pointe, sur le plan de symetrie
Point(3) = {xc, r0, 0};           // fin de l'arrondi, debut de la levre
Point(4) = {xc + Rf, 0, 0};
Point(5) = {xc, Rf, 0};
Point(6) = {xc - Lb, r0, 0};
Point(7) = {xc - Lb, Rf, 0};
// Contour exterieur
Point(10) = {W, 0, 0}; Point(11) = {W, H, 0}; Point(12) = {xf, H, 0}; Point(13) = {xf, r0, 0};
// Trou de goupille
Point(20) = {0, yc, 0}; Point(21) = {R, yc, 0}; Point(22) = {0, yc+R, 0};
Point(23) = {-R, yc, 0}; Point(24) = {0, yc-R, 0};

Circle(1) = {2, 1, 3};   // arrondi de pointe
Circle(2) = {4, 1, 5};   // arc exterieur du bloc A
Line(3) = {2, 4};        // ligament dans le bloc A
Line(4) = {3, 5};        // frontiere A-B
Line(5) = {3, 6};        // levre dans le bloc B
Line(6) = {6, 7};
Line(7) = {5, 7};
Line(8) = {4, 10}; Line(9) = {10, 11}; Line(10) = {11, 12}; Line(11) = {12, 13}; Line(12) = {13, 6};
Circle(20) = {21, 20, 22}; Circle(21) = {22, 20, 23}; Circle(22) = {23, 20, 24}; Circle(23) = {24, 20, 21};

Curve Loop(1) = {3, 2, -4, -1};  Plane Surface(1) = {1};
Curve Loop(2) = {-5, 4, 7, -6};  Plane Surface(2) = {2};
Curve Loop(3) = {8, 9, 10, 11, 12, 6, -7, -2};
Curve Loop(4) = {20, 21, 22, 23};
Plane Surface(3) = {3, 4};

Transfinite Curve {1, 2} = n_theta + 1;
Transfinite Curve {3, 4, 6} = n_r + 1 Using Progression q_r;
Transfinite Curve {5, 7} = n_b + 1 Using Progression q_b;
Transfinite Surface {1} = {2, 4, 5, 3};
Transfinite Surface {2} = {3, 5, 7, 6};
Recombine Surface {1, 2, 3};

Physical Surface("SOLIDE") = {1, 2, 3};
Physical Curve("LIGAMENT") = {3, 8};
Physical Curve("TROU_HAUT") = {20, 21};
Physical Point("BLOC_X") = {10};
Physical Point("POINTE") = {2};

Field[1] = Distance; Field[1].PointsList = {2}; Field[1].Sampling = 100;
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = 0.08; Field[2].SizeMax = h_far;
Field[2].DistMin = 1.5; Field[2].DistMax = 0.4*W;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 1;
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh.MshFileVersion = 2.2;

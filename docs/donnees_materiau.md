# Provenance des données matériau

Toutes les valeurs ont été lues dans l'HDR de B. Tanguy (Université de Lille 1, 2009), qui reproduit les
articles cités. Les numéros de page sont ceux du PDF de l'HDR. Température retenue : **−100 °C**, dans le
bas de la transition ductile-fragile.

## Loi de référence : 16MND5, Tanguy et al. (2005)

R(p, T) = R0 + Q1 (1 − exp(−b1 p)) + Q2 (1 − exp(−b2 p))

| Paramètre | Valeur | Provenance | Fiabilité |
|---|---|---|---|
| E | 210 GPa | Tableau 1 (p. 180) | Valeur publiée |
| ν | 0,3 | Tableau 1 (p. 180) | Valeur publiée |
| Q2 | 472 MPa | Tableau 1 et figure 1 (p. 176, 180) | Valeur publiée |
| b2 | 1,7 | Tableau 1 et figure 1 (p. 176, 180) | Valeur publiée |
| R0(−100 °C) | 589 MPa | Figure 1 (p. 176), lecture en pixels | **Lecture de courbe** |
| Q1(−100 °C) | 184 MPa | Figure 1 (p. 176), lecture en pixels | **Lecture de courbe** |
| b1(−100 °C) | 16,3 | Figure 1 (p. 176), lecture en pixels | **Lecture de courbe** |
| K1(−100 °C) | 84 MPa·s^(1/n1) | Figure 1 (p. 176), lecture en pixels | **Lecture de courbe** |
| n1(−100 °C) | 7,0 | Figure 1 (p. 176), lecture en pixels | **Lecture de courbe** |

### Lecture de la figure 1

La figure a été rendue à 600 dpi. Les axes ont servi d'étalonnage (T de −200 à 200 °C, ordonnée de 0 à 1000 à
gauche et de 0 à 60 à droite), soit environ 1,25 MPa par pixel. Les ordonnées ont été relevées dans la colonne
T = −100 °C. Le marqueur de R0 se trouve à cette température, ce qui peut masquer la courbe de quelques pixels.

**Recoupement :** l'équation (2) de l'article donne une sur-contrainte K1 ε̇^(1/n1) ≈ 22 MPa à 10⁻⁴ s⁻¹, soit
une limite d'élasticité d'environ 611 MPa. La courbe modèle de la figure 2 (p. 176) se lit vers 615 MPa à cette
vitesse, et le point expérimental vers 590 MPa.

### Choix de modélisation

- **Vitesse de déformation :** la loi de référence est indépendante de la vitesse et prend R0 sans
  sur-contrainte. Le cas `tanguy2005_visqueux` ajoute la sur-contrainte de 22 MPa pour en mesurer l'effet.
- **Non pris en compte :** la viscoplasticité complète, l'échauffement adiabatique, l'endommagement ductile.

## Loi de sensibilité : 18MND5, Bordet et al. (2006)

Tableau 5 (p. 116) à −100 °C : E = 208,9 GPa, ν = 0,3, R0 = 539,6 MPa, Q1 = 324,3 MPa, b1 = 12,16,
Q2 = 296,4 MPa, b2 = 0,422, C1/D1 = 73,5 MPa, D1 = 13,5, C2/D2 = 2,2 MPa, D2 = 688,2.

C'est un **autre acier**, utilisé seulement pour mesurer la sensibilité de K_J à la loi d'écrouissage. Ses termes
cinématiques sont cumulés aux termes isotropes, ce qui n'est équivalent qu'en chargement monotone.

## Paramètres de Beremin (pour l'étape 3)

| Source | Acier | Conditions | m | σu (MPa) | V0 |
|---|---|---|---|---|---|
| Tanguy et al. (2005), tableau 1 (p. 180) | 16MND5 | Barreaux entaillés, T ≤ −150 °C | 17,8 | 2925 | 0,001 mm³ |
| Bouchet et al. (2006), figure 28a (p. 96) | 16MND5 | Barreaux entaillés, −100 °C | 19,9 | 3014 | Non trouvé |
| Bouchet et al. (2006), figure 28b (p. 96) | 16MND5 | Barreaux entaillés, −150 à −40 °C | 20,1 | 2976 | Non trouvé |

Choix retenus à l'étape 3 :

- Le jeu de Tanguy et al. est utilisé, car c'est le seul complet et il provient de la même coulée que la loi
  d'écrouissage.
- Ce modèle inclut une correction en déformation plastique : k = 4 et pc = 10⁻⁴ (tableau 1). Une autre page de
  l'HDR donne pc = 0 (p. 200) : la variante pc = 0 change K_J de 0,1 % environ.
- Tanguy et al. (figure 31a, p. 131) gardent σu constant à 2925 MPa jusque vers −80 °C, ce qui justifie cette
  valeur à −100 °C.
- Les isoprobabilités de la figure 31b (p. 131) ont été lues en pixels à −100 °C : 55, 83,5 et 110 MPa√m pour 10,
  50 et 90 % (marge de lecture d'environ 1 MPa√m).
- Le V0 de Bouchet et al. n'est pas indiqué dans les passages lus : **σu ne doit pas être réutilisé sans lui**.
- Ces paramètres dépendent aussi de la loi de comportement, du maillage et de la cinématique utilisés pour les
  identifier. Les réutiliser sur un autre calcul reste une approximation.

## Références

- B. Tanguy, *Fiabilité et durée de vie des composants et des structures métalliques*, HDR, Université de
  Lille 1, 2009. <https://pepite-depot.univ-lille.fr/LIBRE/HDR/2009/50376-2009-379.pdf>
- B. Tanguy, J. Besson, R. Piques, A. Pineau, *Ductile to brittle transition of an A508 steel characterized by
  Charpy impact test. Part II: Modeling of the Charpy transition curve*, Engineering Fracture Mechanics 72
  (2005) 413–434.
- B. Tanguy, C. Bouchet, J. Besson, S. Bugat, *Local approach to fracture based prediction of the ΔT56J and
  ΔTKIc,100 shifts due to irradiation for an A508 pressure vessel steel*, Engineering Fracture Mechanics 73 (2006)
  191–206 (source citée par l'HDR pour la figure 31b).
- A. Pineau, B. Tanguy, *Advances in cleavage fracture modelling in steels*, C. R. Physique 11 (2010) 316–325
  (loi en K⁴ en plasticité confinée, équation 4).
- C. Bouchet, B. Tanguy, J. Besson, A. Pineau, S. Bugat, *Transferability of cleavage fracture parameters between
  notched and cracked geometry*, ECF 16, 2006.
- S. R. Bordet, B. Tanguy, J. Besson, S. Bugat, D. Moinereau, A. Pineau, *Cleavage fracture of RPV steel following
  warm pre-stressing: micromechanical analysis and interpretation through a new model*, Fatigue & Fracture of
  Engineering Materials & Structures, 2006, doi:10.1111/j.1460-2695.2006.01032.x.

## Formules vérifiées

| Formule | Usage ici | État | Source |
|---|---|---|---|
| Charge limite P0 = 1,455 η c σ0, avec η = [(2a/c)² + 2(2a/c) + 2]^½ − (2a/c + 1) | Normaliser la charge | **Confirmée** dans la section sur l'éprouvette CT (éq. 3-6 et 3-8) | V. Kumar, M. D. German, C. F. Shih, rapport EPRI NP-1931, 1981, <https://www.osti.gov/servlets/purl/6068291> |
| Correction d'épaisseur K(x) = 20 + (K(0) − 20)(B0/Bx)^(1/4) | Comparaison de l'étape 4 | **Confirmée** par un article qui cite ASTM E1921 | Li, Chaouadi, Uytdenhouwen, Pardoen, Engineering Fracture Mechanics 290 (2023) 109486 |
| J_pl = η A_pl / (B b) | Contrôle de J par l'aire P–δ | **Forme confirmée** | K. Shi, L. Cai, C. Bao, ICF13, 2013 |
| η_pl = 2 + 0,522 b0/W pour une CT | Contrôle de J par l'aire P–δ | **Coefficient non confirmé** : écrit de mémoire, à vérifier dans une copie légale d'ASTM E1820 | |

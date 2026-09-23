import TzapLean.GF128Proof

/-!
# Rabin certificate for the 128-bit modulus

The residues below are repeated squares of `X` modulo `modulus`. Each step carries an
explicit quotient witness and is checked as an ordinary polynomial identity. The final
residue is `X`; the residue at step 64 also has an explicit Bézout certificate.
-/

namespace TzapLean.GF128Proof

open Polynomial

set_option maxRecDepth 10000
set_option maxHeartbeats 0

@[simp] private theorem coeff2 : (2 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 2
@[simp] private theorem coeff3 : (3 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 3
@[simp] private theorem coeff4 : (4 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 4
@[simp] private theorem coeff5 : (5 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 5
@[simp] private theorem coeff6 : (6 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 6
@[simp] private theorem coeff7 : (7 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 7
@[simp] private theorem coeff8 : (8 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 8
@[simp] private theorem coeff9 : (9 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 9
@[simp] private theorem coeff10 : (10 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 10
@[simp] private theorem coeff11 : (11 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 11
@[simp] private theorem coeff12 : (12 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 12
@[simp] private theorem coeff13 : (13 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 13
@[simp] private theorem coeff14 : (14 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 14
@[simp] private theorem coeff15 : (15 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 15
@[simp] private theorem coeff16 : (16 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 16
@[simp] private theorem coeff17 : (17 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 17
@[simp] private theorem coeff18 : (18 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 18
@[simp] private theorem coeff19 : (19 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 19
@[simp] private theorem coeff20 : (20 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 20
@[simp] private theorem coeff21 : (21 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 21
@[simp] private theorem coeff22 : (22 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 22
@[simp] private theorem coeff23 : (23 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 23
@[simp] private theorem coeff24 : (24 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 24
@[simp] private theorem coeff25 : (25 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 25
@[simp] private theorem coeff26 : (26 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 26
@[simp] private theorem coeff27 : (27 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 27
@[simp] private theorem coeff28 : (28 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 28
@[simp] private theorem coeff29 : (29 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 29
@[simp] private theorem coeff30 : (30 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 30
@[simp] private theorem coeff31 : (31 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 31
@[simp] private theorem coeff32 : (32 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 32
@[simp] private theorem coeff33 : (33 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 33
@[simp] private theorem coeff34 : (34 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 34
@[simp] private theorem coeff35 : (35 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 35
@[simp] private theorem coeff36 : (36 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 36
@[simp] private theorem coeff37 : (37 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 37
@[simp] private theorem coeff38 : (38 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 38
@[simp] private theorem coeff39 : (39 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 39
@[simp] private theorem coeff40 : (40 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 40
@[simp] private theorem coeff41 : (41 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 41
@[simp] private theorem coeff42 : (42 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 42
@[simp] private theorem coeff43 : (43 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 43
@[simp] private theorem coeff44 : (44 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 44
@[simp] private theorem coeff45 : (45 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 45
@[simp] private theorem coeff46 : (46 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 46
@[simp] private theorem coeff47 : (47 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 47
@[simp] private theorem coeff48 : (48 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 48
@[simp] private theorem coeff49 : (49 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 49
@[simp] private theorem coeff50 : (50 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 50
@[simp] private theorem coeff51 : (51 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 51
@[simp] private theorem coeff52 : (52 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 52
@[simp] private theorem coeff53 : (53 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 53
@[simp] private theorem coeff54 : (54 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 54
@[simp] private theorem coeff55 : (55 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 55
@[simp] private theorem coeff56 : (56 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 56
@[simp] private theorem coeff57 : (57 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 57
@[simp] private theorem coeff58 : (58 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 58
@[simp] private theorem coeff59 : (59 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 59
@[simp] private theorem coeff60 : (60 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 60
@[simp] private theorem coeff61 : (61 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 61
@[simp] private theorem coeff62 : (62 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 62
@[simp] private theorem coeff63 : (63 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 63
@[simp] private theorem coeff64 : (64 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 64
@[simp] private theorem coeff65 : (65 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 65
@[simp] private theorem coeff66 : (66 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 66
@[simp] private theorem coeff67 : (67 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 67
@[simp] private theorem coeff68 : (68 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 68
@[simp] private theorem coeff69 : (69 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 69
@[simp] private theorem coeff70 : (70 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 70
@[simp] private theorem coeff71 : (71 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 71
@[simp] private theorem coeff72 : (72 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 72
@[simp] private theorem coeff73 : (73 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 73
@[simp] private theorem coeff74 : (74 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 74
@[simp] private theorem coeff75 : (75 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 75
@[simp] private theorem coeff76 : (76 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 76
@[simp] private theorem coeff77 : (77 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 77
@[simp] private theorem coeff78 : (78 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 78
@[simp] private theorem coeff79 : (79 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 79
@[simp] private theorem coeff80 : (80 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 80
@[simp] private theorem coeff81 : (81 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 81
@[simp] private theorem coeff82 : (82 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 82
@[simp] private theorem coeff83 : (83 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 83
@[simp] private theorem coeff84 : (84 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 84
@[simp] private theorem coeff85 : (85 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 85
@[simp] private theorem coeff86 : (86 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 86
@[simp] private theorem coeff87 : (87 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 87
@[simp] private theorem coeff88 : (88 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 88
@[simp] private theorem coeff89 : (89 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 89
@[simp] private theorem coeff90 : (90 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 90
@[simp] private theorem coeff91 : (91 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 91
@[simp] private theorem coeff92 : (92 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 92
@[simp] private theorem coeff93 : (93 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 93
@[simp] private theorem coeff94 : (94 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 94
@[simp] private theorem coeff95 : (95 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 95
@[simp] private theorem coeff96 : (96 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 96
@[simp] private theorem coeff97 : (97 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 97
@[simp] private theorem coeff98 : (98 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 98
@[simp] private theorem coeff99 : (99 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 99
@[simp] private theorem coeff100 : (100 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 100
@[simp] private theorem coeff101 : (101 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 101
@[simp] private theorem coeff102 : (102 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 102
@[simp] private theorem coeff103 : (103 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 103
@[simp] private theorem coeff104 : (104 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 104
@[simp] private theorem coeff105 : (105 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 105
@[simp] private theorem coeff106 : (106 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 106
@[simp] private theorem coeff107 : (107 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 107
@[simp] private theorem coeff108 : (108 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 108
@[simp] private theorem coeff109 : (109 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 109
@[simp] private theorem coeff110 : (110 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 110
@[simp] private theorem coeff111 : (111 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 111
@[simp] private theorem coeff112 : (112 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 112
@[simp] private theorem coeff113 : (113 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 113
@[simp] private theorem coeff114 : (114 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 114
@[simp] private theorem coeff115 : (115 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 115
@[simp] private theorem coeff116 : (116 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 116
@[simp] private theorem coeff117 : (117 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 117
@[simp] private theorem coeff118 : (118 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 118
@[simp] private theorem coeff119 : (119 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 119
@[simp] private theorem coeff120 : (120 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 120
@[simp] private theorem coeff121 : (121 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 121
@[simp] private theorem coeff122 : (122 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 122
@[simp] private theorem coeff123 : (123 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 123
@[simp] private theorem coeff124 : (124 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 124
@[simp] private theorem coeff125 : (125 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 125
@[simp] private theorem coeff126 : (126 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 126
@[simp] private theorem coeff127 : (127 : F₂[X]) = 1 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 127
@[simp] private theorem coeff128 : (128 : F₂[X]) = 0 := by
  simpa using CharP.cast_eq_mod F₂[X] 2 128

private noncomputable def residue0 : F₂[X] := X
private noncomputable def residue1 : F₂[X] := X ^ 2
private noncomputable def residue2 : F₂[X] := X ^ 4
private noncomputable def residue3 : F₂[X] := X ^ 8
private noncomputable def residue4 : F₂[X] := X ^ 16
private noncomputable def residue5 : F₂[X] := X ^ 32
private noncomputable def residue6 : F₂[X] := X ^ 64
private noncomputable def residue7 : F₂[X] := 1 + X + X ^ 2 + X ^ 7
private noncomputable def residue8 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 14
private noncomputable def residue9 : F₂[X] := 1 + X ^ 4 + X ^ 8 + X ^ 28
private noncomputable def residue10 : F₂[X] := 1 + X ^ 8 + X ^ 16 + X ^ 56
private noncomputable def residue11 : F₂[X] := 1 + X ^ 16 + X ^ 32 + X ^ 112
private noncomputable def residue12 : F₂[X] := 1 + X ^ 32 + X ^ 64 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 103
private noncomputable def residue13 : F₂[X] := X + X ^ 2 + X ^ 7 + X ^ 65 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 73 + X ^ 75 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 85
private noncomputable def residue14 : F₂[X] := X ^ 3 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 39 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 49
private noncomputable def residue15 : F₂[X] := X ^ 6 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 56 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 74 + X ^ 78 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 98
private noncomputable def residue16 : F₂[X] := X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 24 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 35 + X ^ 36 + X ^ 41 + X ^ 42 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 55 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 69 + X ^ 70 + X ^ 72 + X ^ 75 + X ^ 80 + X ^ 84 + X ^ 88 + X ^ 92 + X ^ 96 + X ^ 100 + X ^ 112 + X ^ 124
private noncomputable def residue17 : F₂[X] := 1 + X + X ^ 2 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 19 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 44 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 60 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 70 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 79 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 97 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 110 + X ^ 112 + X ^ 121 + X ^ 122 + X ^ 127
private noncomputable def residue18 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 9 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 18 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 36 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 47 + X ^ 48 + X ^ 53 + X ^ 54 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 72 + X ^ 74 + X ^ 77 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 88 + X ^ 92 + X ^ 93 + X ^ 96 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 103 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 115 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 127
private noncomputable def residue19 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 12 + X ^ 16 + X ^ 17 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 28 + X ^ 29 + X ^ 32 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 52 + X ^ 54 + X ^ 55 + X ^ 57 + X ^ 59 + X ^ 62 + X ^ 63 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 77 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 102 + X ^ 104 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue20 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 29 + X ^ 30 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 119 + X ^ 120
private noncomputable def residue21 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 17 + X ^ 18 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 31 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 43 + X ^ 47 + X ^ 49 + X ^ 55 + X ^ 56 + X ^ 61 + X ^ 65 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 81 + X ^ 84 + X ^ 87 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 104 + X ^ 105 + X ^ 109 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 124
private noncomputable def residue22 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 6 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 39 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 60 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 68 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 127
private noncomputable def residue23 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 39 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 54 + X ^ 55 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 81 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 95 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 108 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 127
private noncomputable def residue24 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 11 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 25 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 44 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 55 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 70 + X ^ 72 + X ^ 77 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 126 + X ^ 127
private noncomputable def residue25 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 11 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 22 + X ^ 23 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 54 + X ^ 56 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 75 + X ^ 76 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 95 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 107 + X ^ 110 + X ^ 114 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue26 : F₂[X] := X ^ 2 + X ^ 6 + X ^ 9 + X ^ 14 + X ^ 15 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 44 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 80 + X ^ 85 + X ^ 86 + X ^ 89 + X ^ 90 + X ^ 96 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 122 + X ^ 123 + X ^ 125
private noncomputable def residue27 : F₂[X] := 1 + X ^ 3 + X ^ 5 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 16 + X ^ 17 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 39 + X ^ 40 + X ^ 43 + X ^ 45 + X ^ 46 + X ^ 49 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 65 + X ^ 68 + X ^ 69 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 101 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 113 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 122 + X ^ 125
private noncomputable def residue28 : F₂[X] := 1 + X + X ^ 4 + X ^ 6 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 61 + X ^ 63 + X ^ 65 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 81 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 97 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 109 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def residue29 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 9 + X ^ 12 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 59 + X ^ 63 + X ^ 65 + X ^ 67 + X ^ 71 + X ^ 73 + X ^ 77 + X ^ 80 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 109 + X ^ 110 + X ^ 113 + X ^ 114 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 126 + X ^ 127
private noncomputable def residue30 : F₂[X] := X + X ^ 6 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 42 + X ^ 45 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 69 + X ^ 70 + X ^ 75 + X ^ 77 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 116 + X ^ 117 + X ^ 120 + X ^ 121 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue31 : F₂[X] := X ^ 4 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 13 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 47 + X ^ 53 + X ^ 54 + X ^ 57 + X ^ 61 + X ^ 62 + X ^ 65 + X ^ 68 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 98 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 125 + X ^ 126
private noncomputable def residue32 : F₂[X] := X + X ^ 3 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 56 + X ^ 61 + X ^ 67 + X ^ 68 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 99 + X ^ 100 + X ^ 103 + X ^ 105 + X ^ 107 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue33 : F₂[X] := X + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 20 + X ^ 24 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 35 + X ^ 37 + X ^ 39 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 49 + X ^ 51 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 80 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 97 + X ^ 101 + X ^ 109 + X ^ 112 + X ^ 113 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 126 + X ^ 127
private noncomputable def residue34 : F₂[X] := X ^ 4 + X ^ 5 + X ^ 6 + X ^ 9 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 49 + X ^ 51 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 73 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 81 + X ^ 84 + X ^ 88 + X ^ 91 + X ^ 96 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 108 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 126 + X ^ 127
private noncomputable def residue35 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 15 + X ^ 20 + X ^ 23 + X ^ 24 + X ^ 30 + X ^ 31 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 69 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 90 + X ^ 91 + X ^ 94 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 108 + X ^ 109 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue36 : F₂[X] := X ^ 2 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 41 + X ^ 43 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 59 + X ^ 67 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 109 + X ^ 116 + X ^ 117 + X ^ 119 + X ^ 125
private noncomputable def residue37 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 6 + X ^ 7 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 45 + X ^ 46 + X ^ 49 + X ^ 51 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 61 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 102 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 113 + X ^ 117 + X ^ 118 + X ^ 122 + X ^ 123 + X ^ 124
private noncomputable def residue38 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 23 + X ^ 24 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 35 + X ^ 37 + X ^ 39 + X ^ 40 + X ^ 42 + X ^ 43 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 52 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 77 + X ^ 85 + X ^ 87 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 117 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue39 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 47 + X ^ 49 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 61 + X ^ 65 + X ^ 69 + X ^ 71 + X ^ 73 + X ^ 76 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 87 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 109 + X ^ 111 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue40 : F₂[X] := X + X ^ 4 + X ^ 9 + X ^ 11 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 21 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 43 + X ^ 46 + X ^ 47 + X ^ 50 + X ^ 52 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 60 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 99 + X ^ 103 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 123 + X ^ 124
private noncomputable def residue41 : F₂[X] := X ^ 3 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 19 + X ^ 22 + X ^ 23 + X ^ 26 + X ^ 27 + X ^ 30 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 49 + X ^ 50 + X ^ 54 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 84 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 95 + X ^ 99 + X ^ 101 + X ^ 102 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 125 + X ^ 127
private noncomputable def residue42 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 7 + X ^ 12 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 27 + X ^ 33 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 78 + X ^ 80 + X ^ 81 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 113 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 126 + X ^ 127
private noncomputable def residue43 : F₂[X] := 1 + X ^ 2 + X ^ 7 + X ^ 10 + X ^ 11 + X ^ 16 + X ^ 17 + X ^ 20 + X ^ 23 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 33 + X ^ 39 + X ^ 43 + X ^ 48 + X ^ 50 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 74 + X ^ 76 + X ^ 80 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue44 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 21 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 31 + X ^ 33 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 67 + X ^ 70 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue45 : F₂[X] := X ^ 3 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 20 + X ^ 21 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 33 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 39 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 65 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 73 + X ^ 74 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 122 + X ^ 125 + X ^ 126
private noncomputable def residue46 : F₂[X] := X + X ^ 3 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 21 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 27 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 75 + X ^ 76 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 91 + X ^ 92 + X ^ 95 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 104 + X ^ 107 + X ^ 113 + X ^ 115 + X ^ 120 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue47 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 9 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 42 + X ^ 43 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 53 + X ^ 57 + X ^ 58 + X ^ 61 + X ^ 62 + X ^ 67 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 77 + X ^ 79 + X ^ 82 + X ^ 86 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 113 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 125 + X ^ 127
private noncomputable def residue48 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 36 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 67 + X ^ 68 + X ^ 73 + X ^ 81 + X ^ 84 + X ^ 87 + X ^ 89 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 105 + X ^ 106 + X ^ 108 + X ^ 109 + X ^ 111 + X ^ 113 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 121 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue49 : F₂[X] := X ^ 6 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 25 + X ^ 30 + X ^ 32 + X ^ 35 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 90 + X ^ 94 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 106 + X ^ 107 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue50 : F₂[X] := X + X ^ 2 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 38 + X ^ 39 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 111 + X ^ 113 + X ^ 116 + X ^ 117 + X ^ 120 + X ^ 122 + X ^ 123 + X ^ 126
private noncomputable def residue51 : F₂[X] := 1 + X + X ^ 3 + X ^ 4 + X ^ 6 + X ^ 7 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 19 + X ^ 25 + X ^ 26 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 49 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 59 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 79 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 99 + X ^ 101 + X ^ 106 + X ^ 107 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 123 + X ^ 126
private noncomputable def residue52 : F₂[X] := X + X ^ 2 + X ^ 4 + X ^ 7 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 22 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 31 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 41 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 60 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 75 + X ^ 77 + X ^ 78 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 97 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 123
private noncomputable def residue53 : F₂[X] := 1 + X + X ^ 5 + X ^ 10 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 37 + X ^ 39 + X ^ 40 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 58 + X ^ 59 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 111 + X ^ 115 + X ^ 118 + X ^ 125 + X ^ 126
private noncomputable def residue54 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 6 + X ^ 12 + X ^ 15 + X ^ 17 + X ^ 21 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 81 + X ^ 84 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 94 + X ^ 95 + X ^ 98 + X ^ 101 + X ^ 103 + X ^ 106 + X ^ 108 + X ^ 115 + X ^ 116 + X ^ 118 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125
private noncomputable def residue55 : F₂[X] := 1 + X + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 13 + X ^ 14 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 27 + X ^ 28 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 70 + X ^ 74 + X ^ 78 + X ^ 79 + X ^ 81 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 95 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 108 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue56 : F₂[X] := 1 + X + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 44 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 52 + X ^ 57 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 68 + X ^ 71 + X ^ 72 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 102 + X ^ 105 + X ^ 108 + X ^ 110 + X ^ 113 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue57 : F₂[X] := X ^ 3 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 36 + X ^ 41 + X ^ 42 + X ^ 47 + X ^ 49 + X ^ 53 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 63 + X ^ 65 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 77 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 100 + X ^ 104 + X ^ 105 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue58 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 47 + X ^ 48 + X ^ 50 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 62 + X ^ 64 + X ^ 69 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 87 + X ^ 89 + X ^ 94 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 114 + X ^ 117 + X ^ 119 + X ^ 122 + X ^ 124
private noncomputable def residue59 : F₂[X] := X + X ^ 2 + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 27 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 38 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 47 + X ^ 51 + X ^ 53 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 97 + X ^ 99 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 127
private noncomputable def residue60 : F₂[X] := X ^ 2 + X ^ 6 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 19 + X ^ 21 + X ^ 23 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 35 + X ^ 36 + X ^ 39 + X ^ 41 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 55 + X ^ 57 + X ^ 58 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 78 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 97 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue61 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 9 + X ^ 11 + X ^ 14 + X ^ 15 + X ^ 19 + X ^ 21 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 29 + X ^ 30 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 73 + X ^ 79 + X ^ 80 + X ^ 83 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 96 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 113 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 127
private noncomputable def residue62 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 6 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 28 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 37 + X ^ 40 + X ^ 42 + X ^ 45 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 54 + X ^ 58 + X ^ 59 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 100 + X ^ 106 + X ^ 107 + X ^ 109 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125
private noncomputable def residue63 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 27 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 36 + X ^ 37 + X ^ 40 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 71 + X ^ 72 + X ^ 79 + X ^ 80 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 101 + X ^ 104 + X ^ 105 + X ^ 107 + X ^ 109 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 118 + X ^ 121 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue64 : F₂[X] := X + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 15 + X ^ 17 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 37 + X ^ 39 + X ^ 43 + X ^ 44 + X ^ 47 + X ^ 48 + X ^ 50 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 97 + X ^ 99 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 125 + X ^ 126
private noncomputable def residue65 : F₂[X] := 1 + X ^ 2 + X ^ 5 + X ^ 6 + X ^ 9 + X ^ 11 + X ^ 13 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 44 + X ^ 47 + X ^ 49 + X ^ 54 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 64 + X ^ 65 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 75 + X ^ 78 + X ^ 79 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 97 + X ^ 98 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 108 + X ^ 109 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 119 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue66 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 13 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 31 + X ^ 32 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 68 + X ^ 70 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 100 + X ^ 101 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 126 + X ^ 127
private noncomputable def residue67 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 7 + X ^ 10 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 21 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 29 + X ^ 30 + X ^ 34 + X ^ 36 + X ^ 41 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 99 + X ^ 101 + X ^ 103 + X ^ 104 + X ^ 106 + X ^ 107 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 117 + X ^ 120 + X ^ 121 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue68 : F₂[X] := X + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 42 + X ^ 45 + X ^ 46 + X ^ 49 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 64 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 84 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 100 + X ^ 103 + X ^ 106 + X ^ 110 + X ^ 115 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue69 : F₂[X] := X + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 41 + X ^ 42 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 78 + X ^ 80 + X ^ 86 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 98 + X ^ 99 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue70 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 33 + X ^ 35 + X ^ 38 + X ^ 39 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 58 + X ^ 62 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 79 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 105 + X ^ 107 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue71 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 67 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 101 + X ^ 102 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 112 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 124 + X ^ 125
private noncomputable def residue72 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 10 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 46 + X ^ 52 + X ^ 56 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 86 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 126 + X ^ 127
private noncomputable def residue73 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 38 + X ^ 41 + X ^ 50 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 57 + X ^ 60 + X ^ 63 + X ^ 66 + X ^ 67 + X ^ 70 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 104 + X ^ 105 + X ^ 113 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 126 + X ^ 127
private noncomputable def residue74 : F₂[X] := 1 + X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 24 + X ^ 26 + X ^ 42 + X ^ 47 + X ^ 54 + X ^ 56 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 87 + X ^ 89 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 106 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue75 : F₂[X] := 1 + X + X ^ 5 + X ^ 9 + X ^ 13 + X ^ 16 + X ^ 18 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 45 + X ^ 46 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 57 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 73 + X ^ 78 + X ^ 79 + X ^ 81 + X ^ 83 + X ^ 85 + X ^ 86 + X ^ 91 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 110 + X ^ 113 + X ^ 116 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 126
private noncomputable def residue76 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 9 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 34 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 54 + X ^ 55 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 81 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 93 + X ^ 94 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue77 : F₂[X] := X + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 61 + X ^ 65 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 74 + X ^ 75 + X ^ 77 + X ^ 81 + X ^ 83 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 94 + X ^ 95 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 105 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 119 + X ^ 122 + X ^ 123 + X ^ 124
private noncomputable def residue78 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 17 + X ^ 19 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 44 + X ^ 47 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 58 + X ^ 60 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 97 + X ^ 99 + X ^ 101 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 125 + X ^ 127
private noncomputable def residue79 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 21 + X ^ 22 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 54 + X ^ 55 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 64 + X ^ 67 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 81 + X ^ 83 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 91 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue80 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 24 + X ^ 31 + X ^ 33 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 47 + X ^ 49 + X ^ 51 + X ^ 53 + X ^ 58 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 72 + X ^ 75 + X ^ 77 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 97 + X ^ 100 + X ^ 101 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue81 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 27 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 51 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 60 + X ^ 64 + X ^ 65 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 72 + X ^ 75 + X ^ 78 + X ^ 79 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 121 + X ^ 123 + X ^ 126 + X ^ 127
private noncomputable def residue82 : F₂[X] := 1 + X ^ 2 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 18 + X ^ 22 + X ^ 26 + X ^ 30 + X ^ 31 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 54 + X ^ 58 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 81 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 99 + X ^ 100 + X ^ 105 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 124 + X ^ 127
private noncomputable def residue83 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 46 + X ^ 47 + X ^ 50 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 83 + X ^ 86 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 101 + X ^ 103 + X ^ 104 + X ^ 106 + X ^ 109 + X ^ 112 + X ^ 114 + X ^ 117 + X ^ 120 + X ^ 122
private noncomputable def residue84 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 27 + X ^ 34 + X ^ 38 + X ^ 43 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 64 + X ^ 67 + X ^ 68 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 87 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 106 + X ^ 110 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 123 + X ^ 124 + X ^ 126
private noncomputable def residue85 : F₂[X] := 1 + X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 21 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 31 + X ^ 34 + X ^ 35 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 47 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 107 + X ^ 111 + X ^ 114 + X ^ 117 + X ^ 118 + X ^ 121 + X ^ 122 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue86 : F₂[X] := 1 + X ^ 4 + X ^ 5 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 52 + X ^ 54 + X ^ 55 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 73 + X ^ 74 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 83 + X ^ 86 + X ^ 87 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 113 + X ^ 114 + X ^ 117 + X ^ 123 + X ^ 125
private noncomputable def residue87 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 5 + X ^ 9 + X ^ 11 + X ^ 13 + X ^ 15 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 40 + X ^ 47 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 62 + X ^ 63 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 93 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 101 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 110 + X ^ 113 + X ^ 114 + X ^ 119 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue88 : F₂[X] := X + X ^ 2 + X ^ 8 + X ^ 9 + X ^ 13 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 37 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 54 + X ^ 56 + X ^ 59 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 77 + X ^ 78 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 110 + X ^ 111 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 124 + X ^ 127
private noncomputable def residue89 : F₂[X] := X ^ 2 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 33 + X ^ 35 + X ^ 36 + X ^ 40 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 53 + X ^ 55 + X ^ 58 + X ^ 59 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 89 + X ^ 91 + X ^ 94 + X ^ 95 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 109 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def residue90 : F₂[X] := 1 + X + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 10 + X ^ 15 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 35 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 47 + X ^ 49 + X ^ 52 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 86 + X ^ 88 + X ^ 91 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue91 : F₂[X] := 1 + X + X ^ 3 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 21 + X ^ 23 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 36 + X ^ 37 + X ^ 41 + X ^ 42 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 67 + X ^ 71 + X ^ 73 + X ^ 75 + X ^ 77 + X ^ 81 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 95 + X ^ 100 + X ^ 106 + X ^ 108 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 118 + X ^ 123 + X ^ 124 + X ^ 125
private noncomputable def residue92 : F₂[X] := 1 + X + X ^ 2 + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 19 + X ^ 21 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 54 + X ^ 59 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 69 + X ^ 73 + X ^ 79 + X ^ 82 + X ^ 85 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 92 + X ^ 94 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue93 : F₂[X] := 1 + X + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 31 + X ^ 32 + X ^ 36 + X ^ 45 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 70 + X ^ 74 + X ^ 77 + X ^ 78 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 90 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 113 + X ^ 118 + X ^ 121 + X ^ 122 + X ^ 127
private noncomputable def residue94 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 6 + X ^ 10 + X ^ 11 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 55 + X ^ 58 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 69 + X ^ 70 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 90 + X ^ 91 + X ^ 96 + X ^ 98 + X ^ 99 + X ^ 105 + X ^ 109 + X ^ 112 + X ^ 114 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 127
private noncomputable def residue95 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 17 + X ^ 19 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 37 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 65 + X ^ 68 + X ^ 69 + X ^ 72 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 83 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 101 + X ^ 103 + X ^ 106 + X ^ 109 + X ^ 112 + X ^ 116 + X ^ 117 + X ^ 122 + X ^ 123 + X ^ 125 + X ^ 126
private noncomputable def residue96 : F₂[X] := 1 + X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 15 + X ^ 21 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 39 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 67 + X ^ 70 + X ^ 71 + X ^ 75 + X ^ 76 + X ^ 79 + X ^ 81 + X ^ 82 + X ^ 88 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 111 + X ^ 113 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 126
private noncomputable def residue97 : F₂[X] := 1 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 15 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 32 + X ^ 34 + X ^ 35 + X ^ 38 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 48 + X ^ 49 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 59 + X ^ 61 + X ^ 64 + X ^ 65 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 71 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 84 + X ^ 88 + X ^ 91 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 105 + X ^ 107 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 125
private noncomputable def residue98 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 19 + X ^ 21 + X ^ 24 + X ^ 25 + X ^ 27 + X ^ 29 + X ^ 30 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 61 + X ^ 65 + X ^ 66 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 81 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 95 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 111 + X ^ 114 + X ^ 115 + X ^ 118 + X ^ 119 + X ^ 123 + X ^ 124
private noncomputable def residue99 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 33 + X ^ 34 + X ^ 37 + X ^ 38 + X ^ 41 + X ^ 42 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 103 + X ^ 107 + X ^ 111 + X ^ 112 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 125 + X ^ 127
private noncomputable def residue100 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 39 + X ^ 40 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 100 + X ^ 101 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 113 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 127
private noncomputable def residue101 : F₂[X] := 1 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 20 + X ^ 21 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 87 + X ^ 88 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 99 + X ^ 102 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 118 + X ^ 121 + X ^ 122 + X ^ 125 + X ^ 127
private noncomputable def residue102 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 6 + X ^ 7 + X ^ 8 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 19 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 63 + X ^ 64 + X ^ 68 + X ^ 69 + X ^ 74 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 109 + X ^ 112 + X ^ 117 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 126 + X ^ 127
private noncomputable def residue103 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 38 + X ^ 42 + X ^ 43 + X ^ 47 + X ^ 52 + X ^ 55 + X ^ 57 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 94 + X ^ 95 + X ^ 100 + X ^ 103 + X ^ 107 + X ^ 108 + X ^ 110 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue104 : F₂[X] := X + X ^ 2 + X ^ 5 + X ^ 9 + X ^ 10 + X ^ 13 + X ^ 23 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 39 + X ^ 45 + X ^ 48 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 95 + X ^ 99 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 123 + X ^ 124
private noncomputable def residue105 : F₂[X] := 1 + X + X ^ 4 + X ^ 6 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 54 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 81 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 95 + X ^ 99 + X ^ 100 + X ^ 106 + X ^ 107 + X ^ 110 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 119 + X ^ 121 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue106 : F₂[X] := 1 + X + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 34 + X ^ 37 + X ^ 38 + X ^ 42 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 61 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue107 : F₂[X] := X ^ 10 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 35 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 44 + X ^ 45 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 71 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 84 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 122 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue108 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 37 + X ^ 38 + X ^ 41 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 71 + X ^ 75 + X ^ 77 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 89 + X ^ 91 + X ^ 94 + X ^ 98 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 110 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue109 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 27 + X ^ 29 + X ^ 35 + X ^ 36 + X ^ 38 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 57 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 93 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 103 + X ^ 104 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 120 + X ^ 122 + X ^ 123
private noncomputable def residue110 : F₂[X] := X + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 37 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 53 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 67 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 81 + X ^ 84 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 95 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 113 + X ^ 115 + X ^ 120 + X ^ 123 + X ^ 124 + X ^ 125
private noncomputable def residue111 : F₂[X] := X + X ^ 3 + X ^ 6 + X ^ 7 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 30 + X ^ 31 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 49 + X ^ 52 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 63 + X ^ 67 + X ^ 72 + X ^ 75 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 82 + X ^ 86 + X ^ 87 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 113 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue112 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 6 + X ^ 7 + X ^ 10 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 35 + X ^ 36 + X ^ 38 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 47 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 57 + X ^ 58 + X ^ 62 + X ^ 65 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 74 + X ^ 75 + X ^ 77 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 104 + X ^ 105 + X ^ 110 + X ^ 113 + X ^ 115 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 126 + X ^ 127
private noncomputable def residue113 : F₂[X] := X + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 24 + X ^ 29 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 43 + X ^ 45 + X ^ 46 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 58 + X ^ 59 + X ^ 60 + X ^ 63 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 84 + X ^ 86 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 93 + X ^ 96 + X ^ 98 + X ^ 103 + X ^ 104 + X ^ 106 + X ^ 107 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 121 + X ^ 123 + X ^ 125 + X ^ 127
private noncomputable def residue114 : F₂[X] := 1 + X ^ 6 + X ^ 7 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 15 + X ^ 17 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 42 + X ^ 44 + X ^ 48 + X ^ 49 + X ^ 54 + X ^ 55 + X ^ 57 + X ^ 60 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 79 + X ^ 80 + X ^ 81 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 97 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 112 + X ^ 113 + X ^ 114 + X ^ 117 + X ^ 119 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 125 + X ^ 127
private noncomputable def residue115 : F₂[X] := X ^ 3 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 28 + X ^ 30 + X ^ 33 + X ^ 34 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 43 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 53 + X ^ 57 + X ^ 64 + X ^ 65 + X ^ 66 + X ^ 69 + X ^ 70 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 90 + X ^ 93 + X ^ 94 + X ^ 98 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 111 + X ^ 112 + X ^ 113 + X ^ 115 + X ^ 119 + X ^ 120 + X ^ 124 + X ^ 125 + X ^ 126
private noncomputable def residue116 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 24 + X ^ 26 + X ^ 32 + X ^ 36 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 58 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 65 + X ^ 66 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 81 + X ^ 82 + X ^ 89 + X ^ 90 + X ^ 91 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 104 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 113 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue117 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 12 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 21 + X ^ 22 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 33 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 42 + X ^ 43 + X ^ 48 + X ^ 50 + X ^ 51 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 67 + X ^ 72 + X ^ 76 + X ^ 77 + X ^ 78 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 89 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 97 + X ^ 101 + X ^ 105 + X ^ 107 + X ^ 108 + X ^ 111 + X ^ 112 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 118 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 127
private noncomputable def residue118 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 17 + X ^ 23 + X ^ 25 + X ^ 27 + X ^ 29 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 42 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 51 + X ^ 53 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 63 + X ^ 64 + X ^ 65 + X ^ 68 + X ^ 69 + X ^ 72 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 87 + X ^ 90 + X ^ 93 + X ^ 94 + X ^ 96 + X ^ 97 + X ^ 98 + X ^ 102 + X ^ 112 + X ^ 113 + X ^ 116 + X ^ 119 + X ^ 124 + X ^ 125
private noncomputable def residue119 : F₂[X] := 1 + X ^ 2 + X ^ 6 + X ^ 7 + X ^ 10 + X ^ 11 + X ^ 15 + X ^ 18 + X ^ 19 + X ^ 20 + X ^ 21 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 45 + X ^ 47 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 71 + X ^ 73 + X ^ 75 + X ^ 77 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 97 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 110 + X ^ 114 + X ^ 117 + X ^ 118 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 126 + X ^ 127
private noncomputable def residue120 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 19 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 33 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 45 + X ^ 47 + X ^ 50 + X ^ 53 + X ^ 57 + X ^ 59 + X ^ 63 + X ^ 64 + X ^ 67 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 72 + X ^ 74 + X ^ 75 + X ^ 80 + X ^ 81 + X ^ 83 + X ^ 85 + X ^ 87 + X ^ 90 + X ^ 92 + X ^ 93 + X ^ 96 + X ^ 99 + X ^ 101 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 109 + X ^ 110 + X ^ 112 + X ^ 117 + X ^ 123
private noncomputable def residue121 : F₂[X] := X + X ^ 2 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 14 + X ^ 16 + X ^ 19 + X ^ 27 + X ^ 29 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 57 + X ^ 58 + X ^ 63 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 75 + X ^ 76 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 93 + X ^ 96 + X ^ 98 + X ^ 99 + X ^ 100 + X ^ 103 + X ^ 107 + X ^ 108 + X ^ 113 + X ^ 114 + X ^ 119 + X ^ 120 + X ^ 125 + X ^ 126
private noncomputable def residue122 : F₂[X] := 1 + X ^ 2 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 14 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 39 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 45 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 53 + X ^ 54 + X ^ 59 + X ^ 60 + X ^ 64 + X ^ 69 + X ^ 70 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 77 + X ^ 78 + X ^ 82 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 105 + X ^ 107 + X ^ 110 + X ^ 111 + X ^ 113 + X ^ 116 + X ^ 117 + X ^ 119 + X ^ 122 + X ^ 123 + X ^ 125
private noncomputable def residue123 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 11 + X ^ 13 + X ^ 16 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 25 + X ^ 28 + X ^ 30 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 50 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 63 + X ^ 65 + X ^ 67 + X ^ 69 + X ^ 70 + X ^ 71 + X ^ 73 + X ^ 75 + X ^ 79 + X ^ 81 + X ^ 84 + X ^ 87 + X ^ 89 + X ^ 90 + X ^ 95 + X ^ 100 + X ^ 101 + X ^ 104 + X ^ 106 + X ^ 107 + X ^ 110 + X ^ 112 + X ^ 113 + X ^ 116 + X ^ 118 + X ^ 119 + X ^ 122 + X ^ 124 + X ^ 125
private noncomputable def residue124 : F₂[X] := X + X ^ 4 + X ^ 7 + X ^ 8 + X ^ 9 + X ^ 10 + X ^ 11 + X ^ 14 + X ^ 15 + X ^ 17 + X ^ 18 + X ^ 20 + X ^ 21 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 26 + X ^ 29 + X ^ 30 + X ^ 31 + X ^ 35 + X ^ 37 + X ^ 48 + X ^ 51 + X ^ 54 + X ^ 56 + X ^ 57 + X ^ 59 + X ^ 60 + X ^ 62 + X ^ 63 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 70 + X ^ 73 + X ^ 74 + X ^ 75 + X ^ 76 + X ^ 79 + X ^ 82 + X ^ 85 + X ^ 91 + X ^ 96 + X ^ 97 + X ^ 102 + X ^ 103 + X ^ 104 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 114 + X ^ 115 + X ^ 116 + X ^ 120 + X ^ 121 + X ^ 122 + X ^ 126 + X ^ 127
private noncomputable def residue125 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 7 + X ^ 9 + X ^ 10 + X ^ 12 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 17 + X ^ 20 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 26 + X ^ 27 + X ^ 28 + X ^ 29 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 61 + X ^ 62 + X ^ 64 + X ^ 65 + X ^ 67 + X ^ 68 + X ^ 70 + X ^ 71 + X ^ 73 + X ^ 74 + X ^ 76 + X ^ 77 + X ^ 79 + X ^ 81 + X ^ 82 + X ^ 83 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 101 + X ^ 102 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 111 + X ^ 113 + X ^ 114 + X ^ 115 + X ^ 117 + X ^ 119 + X ^ 120 + X ^ 121 + X ^ 123 + X ^ 125 + X ^ 126 + X ^ 127
private noncomputable def residue126 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 28 + X ^ 33 + X ^ 35 + X ^ 39 + X ^ 41 + X ^ 45 + X ^ 47 + X ^ 51 + X ^ 53 + X ^ 57 + X ^ 59 + X ^ 63 + X ^ 64 + X ^ 69 + X ^ 70 + X ^ 75 + X ^ 76 + X ^ 81 + X ^ 82 + X ^ 87 + X ^ 88 + X ^ 93 + X ^ 94 + X ^ 97 + X ^ 99 + X ^ 100 + X ^ 103 + X ^ 105 + X ^ 106 + X ^ 109 + X ^ 111 + X ^ 112 + X ^ 115 + X ^ 117 + X ^ 118 + X ^ 121 + X ^ 123 + X ^ 124 + X ^ 127
private noncomputable def residue127 : F₂[X] := X ^ 2 + X ^ 5 + X ^ 7 + X ^ 8 + X ^ 10 + X ^ 11 + X ^ 13 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 19 + X ^ 20 + X ^ 22 + X ^ 23 + X ^ 25 + X ^ 26 + X ^ 28 + X ^ 29 + X ^ 31 + X ^ 32 + X ^ 34 + X ^ 35 + X ^ 37 + X ^ 38 + X ^ 40 + X ^ 41 + X ^ 43 + X ^ 44 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 52 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 58 + X ^ 59 + X ^ 61 + X ^ 62 + X ^ 65 + X ^ 68 + X ^ 71 + X ^ 74 + X ^ 77 + X ^ 80 + X ^ 83 + X ^ 86 + X ^ 89 + X ^ 92 + X ^ 95 + X ^ 98 + X ^ 101 + X ^ 104 + X ^ 107 + X ^ 110 + X ^ 113 + X ^ 116 + X ^ 119 + X ^ 122 + X ^ 125
private noncomputable def residue128 : F₂[X] := X

private noncomputable def quotient0 : F₂[X] := 0
private noncomputable def quotient1 : F₂[X] := 0
private noncomputable def quotient2 : F₂[X] := 0
private noncomputable def quotient3 : F₂[X] := 0
private noncomputable def quotient4 : F₂[X] := 0
private noncomputable def quotient5 : F₂[X] := 0
private noncomputable def quotient6 : F₂[X] := 1
private noncomputable def quotient7 : F₂[X] := 0
private noncomputable def quotient8 : F₂[X] := 0
private noncomputable def quotient9 : F₂[X] := 0
private noncomputable def quotient10 : F₂[X] := 0
private noncomputable def quotient11 : F₂[X] := X ^ 96
private noncomputable def quotient12 : F₂[X] := 1 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 78
private noncomputable def quotient13 : F₂[X] := X ^ 2 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 18 + X ^ 22 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 42
private noncomputable def quotient14 : F₂[X] := 0
private noncomputable def quotient15 : F₂[X] := X ^ 4 + X ^ 8 + X ^ 12 + X ^ 20 + X ^ 28 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 68
private noncomputable def quotient16 : F₂[X] := 1 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 22 + X ^ 32 + X ^ 40 + X ^ 48 + X ^ 56 + X ^ 64 + X ^ 72 + X ^ 96 + X ^ 120
private noncomputable def quotient17 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 12 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 30 + X ^ 36 + X ^ 40 + X ^ 52 + X ^ 56 + X ^ 60 + X ^ 66 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 92 + X ^ 96 + X ^ 114 + X ^ 116 + X ^ 126
private noncomputable def quotient18 : F₂[X] := X ^ 5 + X ^ 16 + X ^ 20 + X ^ 26 + X ^ 28 + X ^ 32 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 48 + X ^ 56 + X ^ 58 + X ^ 64 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 78 + X ^ 88 + X ^ 92 + X ^ 96 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 126
private noncomputable def quotient19 : F₂[X] := 1 + X + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 16 + X ^ 20 + X ^ 26 + X ^ 28 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 76 + X ^ 80 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 96 + X ^ 100 + X ^ 104 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient20 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 56 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 110 + X ^ 112
private noncomputable def quotient21 : F₂[X] := X ^ 2 + X ^ 12 + X ^ 16 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 32 + X ^ 34 + X ^ 40 + X ^ 46 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 80 + X ^ 82 + X ^ 90 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 120
private noncomputable def quotient22 : F₂[X] := X ^ 2 + X ^ 5 + X ^ 8 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 88 + X ^ 96 + X ^ 100 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 126
private noncomputable def quotient23 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 34 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 56 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 88 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 126
private noncomputable def quotient24 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 12 + X ^ 16 + X ^ 26 + X ^ 28 + X ^ 32 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 124 + X ^ 126
private noncomputable def quotient25 : F₂[X] := 1 + X + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 22 + X ^ 24 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 48 + X ^ 52 + X ^ 56 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 86 + X ^ 92 + X ^ 100 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient26 : F₂[X] := 1 + X + X ^ 4 + X ^ 6 + X ^ 16 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 32 + X ^ 42 + X ^ 44 + X ^ 50 + X ^ 52 + X ^ 64 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 116 + X ^ 118 + X ^ 122
private noncomputable def quotient27 : F₂[X] := X + X ^ 2 + X ^ 8 + X ^ 10 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 74 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 116 + X ^ 122
private noncomputable def quotient28 : F₂[X] := X ^ 2 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 32 + X ^ 34 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 64 + X ^ 66 + X ^ 76 + X ^ 84 + X ^ 88 + X ^ 90 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 112 + X ^ 116 + X ^ 120
private noncomputable def quotient29 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 14 + X ^ 18 + X ^ 26 + X ^ 32 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 90 + X ^ 92 + X ^ 98 + X ^ 100 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 124 + X ^ 126
private noncomputable def quotient30 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 12 + X ^ 22 + X ^ 26 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 104 + X ^ 106 + X ^ 112 + X ^ 114 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient31 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 8 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 56 + X ^ 60 + X ^ 62 + X ^ 68 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 122 + X ^ 124
private noncomputable def quotient32 : F₂[X] := X + X ^ 3 + X ^ 6 + X ^ 8 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 78 + X ^ 82 + X ^ 86 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient33 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 32 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 66 + X ^ 74 + X ^ 90 + X ^ 96 + X ^ 98 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 124 + X ^ 126
private noncomputable def quotient34 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 18 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 34 + X ^ 40 + X ^ 48 + X ^ 54 + X ^ 64 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 88 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 124 + X ^ 126
private noncomputable def quotient35 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 52 + X ^ 54 + X ^ 60 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 88 + X ^ 90 + X ^ 96 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient36 : F₂[X] := X + X ^ 6 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 90 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 122
private noncomputable def quotient37 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 6 + X ^ 10 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 76 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 98 + X ^ 106 + X ^ 108 + X ^ 116 + X ^ 118 + X ^ 120
private noncomputable def quotient38 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 42 + X ^ 46 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 106 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient39 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 10 + X ^ 14 + X ^ 18 + X ^ 24 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 46 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 90 + X ^ 94 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient40 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 6 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 70 + X ^ 78 + X ^ 84 + X ^ 88 + X ^ 92 + X ^ 100 + X ^ 104 + X ^ 108 + X ^ 118 + X ^ 120
private noncomputable def quotient41 : F₂[X] := 1 + X + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 40 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 62 + X ^ 70 + X ^ 74 + X ^ 76 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 122 + X ^ 126
private noncomputable def quotient42 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 28 + X ^ 32 + X ^ 34 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 98 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 124 + X ^ 126
private noncomputable def quotient43 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 20 + X ^ 24 + X ^ 32 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient44 : F₂[X] := X ^ 3 + X ^ 5 + X ^ 6 + X ^ 12 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient45 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 18 + X ^ 20 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 80 + X ^ 84 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 116 + X ^ 122 + X ^ 124
private noncomputable def quotient46 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 12 + X ^ 16 + X ^ 22 + X ^ 24 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 54 + X ^ 56 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 80 + X ^ 86 + X ^ 98 + X ^ 102 + X ^ 112 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient47 : F₂[X] := 1 + X + X ^ 5 + X ^ 6 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 26 + X ^ 30 + X ^ 36 + X ^ 44 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 98 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 122 + X ^ 126
private noncomputable def quotient48 : F₂[X] := 1 + X + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 18 + X ^ 34 + X ^ 40 + X ^ 46 + X ^ 50 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 82 + X ^ 84 + X ^ 88 + X ^ 90 + X ^ 94 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 112 + X ^ 114 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient49 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 30 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 52 + X ^ 60 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 84 + X ^ 86 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient50 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 94 + X ^ 98 + X ^ 104 + X ^ 106 + X ^ 112 + X ^ 116 + X ^ 118 + X ^ 124
private noncomputable def quotient51 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 30 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 70 + X ^ 74 + X ^ 84 + X ^ 86 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 116 + X ^ 118 + X ^ 124
private noncomputable def quotient52 : F₂[X] := 1 + X ^ 4 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 22 + X ^ 26 + X ^ 28 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 118
private noncomputable def quotient53 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 10 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 72 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 94 + X ^ 102 + X ^ 108 + X ^ 122 + X ^ 124
private noncomputable def quotient54 : F₂[X] := X + X ^ 2 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 34 + X ^ 40 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 60 + X ^ 62 + X ^ 68 + X ^ 74 + X ^ 78 + X ^ 84 + X ^ 88 + X ^ 102 + X ^ 104 + X ^ 108 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122
private noncomputable def quotient55 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 12 + X ^ 20 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 62 + X ^ 64 + X ^ 68 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 88 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient56 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 8 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 52 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 76 + X ^ 82 + X ^ 88 + X ^ 92 + X ^ 98 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient57 : F₂[X] := 1 + X + X ^ 2 + X ^ 5 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 26 + X ^ 28 + X ^ 36 + X ^ 40 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 72 + X ^ 80 + X ^ 82 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient58 : F₂[X] := 1 + X ^ 10 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 46 + X ^ 50 + X ^ 60 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 100 + X ^ 106 + X ^ 110 + X ^ 116 + X ^ 120
private noncomputable def quotient59 : F₂[X] := X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 70 + X ^ 76 + X ^ 84 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 126
private noncomputable def quotient60 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 28 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 66 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient61 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 12 + X ^ 16 + X ^ 18 + X ^ 30 + X ^ 32 + X ^ 38 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 64 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 126
private noncomputable def quotient62 : F₂[X] := 1 + X + X ^ 4 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 72 + X ^ 84 + X ^ 86 + X ^ 90 + X ^ 92 + X ^ 100 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122
private noncomputable def quotient63 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 14 + X ^ 16 + X ^ 30 + X ^ 32 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 74 + X ^ 80 + X ^ 82 + X ^ 86 + X ^ 90 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 108 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient64 : F₂[X] := 1 + X + X ^ 3 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 66 + X ^ 70 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 96 + X ^ 100 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 122 + X ^ 124
private noncomputable def quotient65 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 22 + X ^ 28 + X ^ 30 + X ^ 44 + X ^ 48 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 88 + X ^ 90 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 110 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient66 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 8 + X ^ 12 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 72 + X ^ 74 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 100 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 124 + X ^ 126
private noncomputable def quotient67 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 52 + X ^ 56 + X ^ 60 + X ^ 62 + X ^ 70 + X ^ 74 + X ^ 78 + X ^ 80 + X ^ 84 + X ^ 86 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 106 + X ^ 112 + X ^ 114 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient68 : F₂[X] := X + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 40 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 72 + X ^ 78 + X ^ 84 + X ^ 92 + X ^ 102 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient69 : F₂[X] := 1 + X + X ^ 3 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 28 + X ^ 32 + X ^ 44 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 68 + X ^ 70 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 108 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient70 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 30 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 72 + X ^ 82 + X ^ 86 + X ^ 92 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient71 : F₂[X] := 1 + X + X ^ 6 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 60 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 74 + X ^ 76 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 120 + X ^ 122
private noncomputable def quotient72 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 44 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 124 + X ^ 126
private noncomputable def quotient73 : F₂[X] := 1 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 12 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 80 + X ^ 82 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 124 + X ^ 126
private noncomputable def quotient74 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 46 + X ^ 50 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 84 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient75 : F₂[X] := X ^ 3 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 18 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 38 + X ^ 42 + X ^ 44 + X ^ 54 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 92 + X ^ 98 + X ^ 104 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 124
private noncomputable def quotient76 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 12 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 36 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 58 + X ^ 60 + X ^ 68 + X ^ 76 + X ^ 80 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient77 : F₂[X] := X ^ 2 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 34 + X ^ 38 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 60 + X ^ 62 + X ^ 68 + X ^ 72 + X ^ 80 + X ^ 82 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 110 + X ^ 116 + X ^ 118 + X ^ 120
private noncomputable def quotient78 : F₂[X] := X + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 66 + X ^ 70 + X ^ 74 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 122 + X ^ 126
private noncomputable def quotient79 : F₂[X] := X + X ^ 5 + X ^ 6 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 34 + X ^ 38 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 54 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient80 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 22 + X ^ 26 + X ^ 30 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 66 + X ^ 72 + X ^ 74 + X ^ 80 + X ^ 84 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient81 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 16 + X ^ 22 + X ^ 28 + X ^ 30 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 114 + X ^ 118 + X ^ 124 + X ^ 126
private noncomputable def quotient82 : F₂[X] := X ^ 2 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 16 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 34 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 56 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 82 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 120 + X ^ 126
private noncomputable def quotient83 : F₂[X] := 1 + X ^ 2 + X ^ 4 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 38 + X ^ 44 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 72 + X ^ 74 + X ^ 78 + X ^ 80 + X ^ 84 + X ^ 90 + X ^ 96 + X ^ 100 + X ^ 106 + X ^ 112 + X ^ 116
private noncomputable def quotient84 : F₂[X] := 1 + X ^ 3 + X ^ 6 + X ^ 8 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 46 + X ^ 52 + X ^ 60 + X ^ 64 + X ^ 68 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 84 + X ^ 92 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 118 + X ^ 120 + X ^ 124
private noncomputable def quotient85 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 86 + X ^ 94 + X ^ 100 + X ^ 106 + X ^ 108 + X ^ 114 + X ^ 116 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient86 : F₂[X] := 1 + X + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 18 + X ^ 20 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 38 + X ^ 44 + X ^ 46 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 72 + X ^ 76 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 98 + X ^ 100 + X ^ 106 + X ^ 118 + X ^ 122
private noncomputable def quotient87 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 58 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 74 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 92 + X ^ 98 + X ^ 100 + X ^ 110 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient88 : F₂[X] := X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 26 + X ^ 28 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 92 + X ^ 94 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 120 + X ^ 126
private noncomputable def quotient89 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 8 + X ^ 10 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 50 + X ^ 54 + X ^ 60 + X ^ 62 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 90 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 116 + X ^ 120 + X ^ 124
private noncomputable def quotient90 : F₂[X] := X + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 44 + X ^ 48 + X ^ 54 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 76 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient91 : F₂[X] := X + X ^ 2 + X ^ 6 + X ^ 14 + X ^ 18 + X ^ 22 + X ^ 26 + X ^ 34 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 62 + X ^ 72 + X ^ 84 + X ^ 88 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 118 + X ^ 120 + X ^ 122
private noncomputable def quotient92 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 18 + X ^ 30 + X ^ 36 + X ^ 42 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 56 + X ^ 60 + X ^ 76 + X ^ 80 + X ^ 84 + X ^ 92 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient93 : F₂[X] := X ^ 4 + X ^ 5 + X ^ 6 + X ^ 12 + X ^ 20 + X ^ 26 + X ^ 28 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 52 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 76 + X ^ 84 + X ^ 98 + X ^ 108 + X ^ 114 + X ^ 116 + X ^ 126
private noncomputable def quotient94 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 5 + X ^ 10 + X ^ 12 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 52 + X ^ 54 + X ^ 64 + X ^ 68 + X ^ 70 + X ^ 82 + X ^ 90 + X ^ 96 + X ^ 100 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 126
private noncomputable def quotient95 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 8 + X ^ 10 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 38 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 74 + X ^ 78 + X ^ 84 + X ^ 90 + X ^ 96 + X ^ 104 + X ^ 106 + X ^ 116 + X ^ 118 + X ^ 122 + X ^ 124
private noncomputable def quotient96 : F₂[X] := X ^ 2 + X ^ 3 + X ^ 6 + X ^ 12 + X ^ 14 + X ^ 22 + X ^ 24 + X ^ 30 + X ^ 34 + X ^ 36 + X ^ 48 + X ^ 52 + X ^ 60 + X ^ 64 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 94 + X ^ 98 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 124
private noncomputable def quotient97 : F₂[X] := 1 + X + X ^ 2 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 40 + X ^ 48 + X ^ 54 + X ^ 64 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 82 + X ^ 86 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 122
private noncomputable def quotient98 : F₂[X] := X ^ 2 + X ^ 4 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 62 + X ^ 64 + X ^ 68 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 94 + X ^ 100 + X ^ 102 + X ^ 108 + X ^ 110 + X ^ 118 + X ^ 120
private noncomputable def quotient99 : F₂[X] := 1 + X + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 64 + X ^ 68 + X ^ 72 + X ^ 78 + X ^ 86 + X ^ 94 + X ^ 96 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 122 + X ^ 126
private noncomputable def quotient100 : F₂[X] := X ^ 5 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 18 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 72 + X ^ 74 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 98 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 126
private noncomputable def quotient101 : F₂[X] := 1 + X + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 46 + X ^ 48 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 70 + X ^ 76 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 108 + X ^ 114 + X ^ 116 + X ^ 122 + X ^ 126
private noncomputable def quotient102 : F₂[X] := X ^ 3 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 20 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 72 + X ^ 76 + X ^ 84 + X ^ 88 + X ^ 90 + X ^ 96 + X ^ 106 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 124 + X ^ 126
private noncomputable def quotient103 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 60 + X ^ 62 + X ^ 72 + X ^ 78 + X ^ 86 + X ^ 88 + X ^ 92 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient104 : F₂[X] := 1 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 28 + X ^ 32 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 62 + X ^ 70 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 118 + X ^ 120
private noncomputable def quotient105 : F₂[X] := X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 62 + X ^ 70 + X ^ 72 + X ^ 84 + X ^ 86 + X ^ 92 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 110 + X ^ 114 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient106 : F₂[X] := 1 + X + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 64 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient107 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 6 + X ^ 8 + X ^ 14 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 40 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 68 + X ^ 76 + X ^ 80 + X ^ 84 + X ^ 88 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 116 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient108 : F₂[X] := X + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 14 + X ^ 22 + X ^ 26 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 50 + X ^ 54 + X ^ 60 + X ^ 68 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 92 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient109 : F₂[X] := 1 + X ^ 4 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 58 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 78 + X ^ 80 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 100 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 112 + X ^ 116 + X ^ 118
private noncomputable def quotient110 : F₂[X] := X + X ^ 6 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 34 + X ^ 40 + X ^ 46 + X ^ 48 + X ^ 52 + X ^ 62 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 76 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 98 + X ^ 102 + X ^ 112 + X ^ 118 + X ^ 120 + X ^ 122
private noncomputable def quotient111 : F₂[X] := 1 + X + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 16 + X ^ 22 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 36 + X ^ 44 + X ^ 46 + X ^ 54 + X ^ 56 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 98 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient112 : F₂[X] := 1 + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 36 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 56 + X ^ 64 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 80 + X ^ 82 + X ^ 92 + X ^ 98 + X ^ 102 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 124 + X ^ 126
private noncomputable def quotient113 : F₂[X] := 1 + X + X ^ 2 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 16 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 40 + X ^ 44 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 58 + X ^ 64 + X ^ 68 + X ^ 78 + X ^ 80 + X ^ 84 + X ^ 86 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 104 + X ^ 106 + X ^ 108 + X ^ 110 + X ^ 114 + X ^ 118 + X ^ 122 + X ^ 126
private noncomputable def quotient114 : F₂[X] := 1 + X + X ^ 5 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 28 + X ^ 30 + X ^ 32 + X ^ 34 + X ^ 36 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 52 + X ^ 54 + X ^ 56 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 98 + X ^ 100 + X ^ 106 + X ^ 110 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 122 + X ^ 126
private noncomputable def quotient115 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 10 + X ^ 12 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 32 + X ^ 36 + X ^ 40 + X ^ 44 + X ^ 48 + X ^ 50 + X ^ 52 + X ^ 58 + X ^ 60 + X ^ 68 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 94 + X ^ 96 + X ^ 98 + X ^ 102 + X ^ 110 + X ^ 112 + X ^ 120 + X ^ 122 + X ^ 124
private noncomputable def quotient116 : F₂[X] := 1 + X + X ^ 2 + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 34 + X ^ 36 + X ^ 50 + X ^ 52 + X ^ 54 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 80 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 94 + X ^ 98 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient117 : F₂[X] := X ^ 2 + X ^ 5 + X ^ 6 + X ^ 16 + X ^ 24 + X ^ 26 + X ^ 28 + X ^ 30 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 50 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 66 + X ^ 74 + X ^ 82 + X ^ 86 + X ^ 88 + X ^ 94 + X ^ 96 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 126
private noncomputable def quotient118 : F₂[X] := 1 + X + X ^ 2 + X ^ 8 + X ^ 10 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 46 + X ^ 52 + X ^ 58 + X ^ 60 + X ^ 64 + X ^ 66 + X ^ 68 + X ^ 76 + X ^ 96 + X ^ 98 + X ^ 104 + X ^ 110 + X ^ 120 + X ^ 122
private noncomputable def quotient119 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 14 + X ^ 18 + X ^ 22 + X ^ 26 + X ^ 36 + X ^ 38 + X ^ 40 + X ^ 52 + X ^ 56 + X ^ 66 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 92 + X ^ 100 + X ^ 106 + X ^ 108 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 118 + X ^ 120 + X ^ 124 + X ^ 126
private noncomputable def quotient120 : F₂[X] := 1 + X ^ 6 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 20 + X ^ 22 + X ^ 32 + X ^ 34 + X ^ 38 + X ^ 42 + X ^ 46 + X ^ 52 + X ^ 56 + X ^ 58 + X ^ 64 + X ^ 70 + X ^ 74 + X ^ 76 + X ^ 80 + X ^ 84 + X ^ 90 + X ^ 92 + X ^ 96 + X ^ 106 + X ^ 118
private noncomputable def quotient121 : F₂[X] := 1 + X + X ^ 3 + X ^ 12 + X ^ 16 + X ^ 22 + X ^ 24 + X ^ 38 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 58 + X ^ 64 + X ^ 68 + X ^ 70 + X ^ 72 + X ^ 78 + X ^ 86 + X ^ 88 + X ^ 98 + X ^ 100 + X ^ 110 + X ^ 112 + X ^ 122 + X ^ 124
private noncomputable def quotient122 : F₂[X] := 1 + X + X ^ 10 + X ^ 12 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 26 + X ^ 28 + X ^ 36 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 82 + X ^ 86 + X ^ 92 + X ^ 94 + X ^ 98 + X ^ 104 + X ^ 106 + X ^ 110 + X ^ 116 + X ^ 118 + X ^ 122
private noncomputable def quotient123 : F₂[X] := X + X ^ 2 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 18 + X ^ 22 + X ^ 30 + X ^ 34 + X ^ 40 + X ^ 46 + X ^ 50 + X ^ 52 + X ^ 62 + X ^ 72 + X ^ 74 + X ^ 80 + X ^ 84 + X ^ 86 + X ^ 92 + X ^ 96 + X ^ 98 + X ^ 104 + X ^ 108 + X ^ 110 + X ^ 116 + X ^ 120 + X ^ 122
private noncomputable def quotient124 : F₂[X] := X ^ 3 + X ^ 4 + X ^ 5 + X ^ 8 + X ^ 10 + X ^ 12 + X ^ 18 + X ^ 20 + X ^ 22 + X ^ 24 + X ^ 30 + X ^ 36 + X ^ 42 + X ^ 54 + X ^ 64 + X ^ 66 + X ^ 76 + X ^ 78 + X ^ 80 + X ^ 88 + X ^ 90 + X ^ 92 + X ^ 100 + X ^ 102 + X ^ 104 + X ^ 112 + X ^ 114 + X ^ 116 + X ^ 124 + X ^ 126
private noncomputable def quotient125 : F₂[X] := X + X ^ 2 + X ^ 3 + X ^ 5 + X ^ 6 + X ^ 8 + X ^ 12 + X ^ 14 + X ^ 18 + X ^ 20 + X ^ 24 + X ^ 26 + X ^ 30 + X ^ 34 + X ^ 36 + X ^ 38 + X ^ 42 + X ^ 46 + X ^ 48 + X ^ 50 + X ^ 54 + X ^ 58 + X ^ 60 + X ^ 62 + X ^ 64 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 74 + X ^ 76 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 86 + X ^ 88 + X ^ 90 + X ^ 94 + X ^ 98 + X ^ 100 + X ^ 102 + X ^ 106 + X ^ 110 + X ^ 112 + X ^ 114 + X ^ 118 + X ^ 122 + X ^ 124 + X ^ 126
private noncomputable def quotient126 : F₂[X] := X ^ 5 + X ^ 10 + X ^ 12 + X ^ 22 + X ^ 24 + X ^ 34 + X ^ 36 + X ^ 46 + X ^ 48 + X ^ 58 + X ^ 60 + X ^ 66 + X ^ 70 + X ^ 72 + X ^ 78 + X ^ 82 + X ^ 84 + X ^ 90 + X ^ 94 + X ^ 96 + X ^ 102 + X ^ 106 + X ^ 108 + X ^ 114 + X ^ 118 + X ^ 120 + X ^ 126
private noncomputable def quotient127 : F₂[X] := X + X ^ 2 + X ^ 8 + X ^ 14 + X ^ 20 + X ^ 26 + X ^ 32 + X ^ 38 + X ^ 44 + X ^ 50 + X ^ 56 + X ^ 62 + X ^ 68 + X ^ 74 + X ^ 80 + X ^ 86 + X ^ 92 + X ^ 98 + X ^ 104 + X ^ 110 + X ^ 116 + X ^ 122

private theorem square0 : residue0 ^ 2 = residue1 + modulus * quotient0 := by
  unfold residue0 residue1 quotient0 modulus
  ring_nf
private theorem square1 : residue1 ^ 2 = residue2 + modulus * quotient1 := by
  unfold residue1 residue2 quotient1 modulus
  ring_nf
private theorem square2 : residue2 ^ 2 = residue3 + modulus * quotient2 := by
  unfold residue2 residue3 quotient2 modulus
  ring_nf
private theorem square3 : residue3 ^ 2 = residue4 + modulus * quotient3 := by
  unfold residue3 residue4 quotient3 modulus
  ring_nf
private theorem square4 : residue4 ^ 2 = residue5 + modulus * quotient4 := by
  unfold residue4 residue5 quotient4 modulus
  ring_nf
private theorem square5 : residue5 ^ 2 = residue6 + modulus * quotient5 := by
  unfold residue5 residue6 quotient5 modulus
  ring_nf
private theorem square6 : residue6 ^ 2 = residue7 + modulus * quotient6 := by
  unfold residue6 residue7 quotient6 modulus
  ring_nf
  all_goals simp
private theorem square7 : residue7 ^ 2 = residue8 + modulus * quotient7 := by
  unfold residue7 residue8 quotient7 modulus
  ring_nf
  all_goals simp
private theorem square8 : residue8 ^ 2 = residue9 + modulus * quotient8 := by
  unfold residue8 residue9 quotient8 modulus
  ring_nf
  all_goals simp
private theorem square9 : residue9 ^ 2 = residue10 + modulus * quotient9 := by
  unfold residue9 residue10 quotient9 modulus
  ring_nf
  all_goals simp
private theorem square10 : residue10 ^ 2 = residue11 + modulus * quotient10 := by
  unfold residue10 residue11 quotient10 modulus
  ring_nf
  all_goals simp
private theorem square11 : residue11 ^ 2 = residue12 + modulus * quotient11 := by
  unfold residue11 residue12 quotient11 modulus
  ring_nf
  all_goals simp
private theorem square12 : residue12 ^ 2 = residue13 + modulus * quotient12 := by
  unfold residue12 residue13 quotient12 modulus
  ring_nf
  all_goals simp
private theorem square13 : residue13 ^ 2 = residue14 + modulus * quotient13 := by
  unfold residue13 residue14 quotient13 modulus
  ring_nf
  all_goals simp
private theorem square14 : residue14 ^ 2 = residue15 + modulus * quotient14 := by
  unfold residue14 residue15 quotient14 modulus
  ring_nf
  all_goals simp
private theorem square15 : residue15 ^ 2 = residue16 + modulus * quotient15 := by
  unfold residue15 residue16 quotient15 modulus
  ring_nf
  all_goals simp
private theorem square16 : residue16 ^ 2 = residue17 + modulus * quotient16 := by
  unfold residue16 residue17 quotient16 modulus
  ring_nf
  all_goals simp
private theorem square17 : residue17 ^ 2 = residue18 + modulus * quotient17 := by
  unfold residue17 residue18 quotient17 modulus
  ring_nf
  all_goals simp
private theorem square18 : residue18 ^ 2 = residue19 + modulus * quotient18 := by
  unfold residue18 residue19 quotient18 modulus
  ring_nf
  all_goals simp
private theorem square19 : residue19 ^ 2 = residue20 + modulus * quotient19 := by
  unfold residue19 residue20 quotient19 modulus
  ring_nf
  all_goals simp
private theorem square20 : residue20 ^ 2 = residue21 + modulus * quotient20 := by
  unfold residue20 residue21 quotient20 modulus
  ring_nf
  all_goals simp
private theorem square21 : residue21 ^ 2 = residue22 + modulus * quotient21 := by
  unfold residue21 residue22 quotient21 modulus
  ring_nf
  all_goals simp
private theorem square22 : residue22 ^ 2 = residue23 + modulus * quotient22 := by
  unfold residue22 residue23 quotient22 modulus
  ring_nf
  all_goals simp
private theorem square23 : residue23 ^ 2 = residue24 + modulus * quotient23 := by
  unfold residue23 residue24 quotient23 modulus
  ring_nf
  all_goals simp
private theorem square24 : residue24 ^ 2 = residue25 + modulus * quotient24 := by
  unfold residue24 residue25 quotient24 modulus
  ring_nf
  all_goals simp
private theorem square25 : residue25 ^ 2 = residue26 + modulus * quotient25 := by
  unfold residue25 residue26 quotient25 modulus
  ring_nf
  all_goals simp
private theorem square26 : residue26 ^ 2 = residue27 + modulus * quotient26 := by
  unfold residue26 residue27 quotient26 modulus
  ring_nf
  all_goals simp
private theorem square27 : residue27 ^ 2 = residue28 + modulus * quotient27 := by
  unfold residue27 residue28 quotient27 modulus
  ring_nf
  all_goals simp
private theorem square28 : residue28 ^ 2 = residue29 + modulus * quotient28 := by
  unfold residue28 residue29 quotient28 modulus
  ring_nf
  all_goals simp
private theorem square29 : residue29 ^ 2 = residue30 + modulus * quotient29 := by
  unfold residue29 residue30 quotient29 modulus
  ring_nf
  all_goals simp
private theorem square30 : residue30 ^ 2 = residue31 + modulus * quotient30 := by
  unfold residue30 residue31 quotient30 modulus
  ring_nf
  all_goals simp
private theorem square31 : residue31 ^ 2 = residue32 + modulus * quotient31 := by
  unfold residue31 residue32 quotient31 modulus
  ring_nf
  all_goals simp
private theorem square32 : residue32 ^ 2 = residue33 + modulus * quotient32 := by
  unfold residue32 residue33 quotient32 modulus
  ring_nf
  all_goals simp
private theorem square33 : residue33 ^ 2 = residue34 + modulus * quotient33 := by
  unfold residue33 residue34 quotient33 modulus
  ring_nf
  all_goals simp
private theorem square34 : residue34 ^ 2 = residue35 + modulus * quotient34 := by
  unfold residue34 residue35 quotient34 modulus
  ring_nf
  all_goals simp
private theorem square35 : residue35 ^ 2 = residue36 + modulus * quotient35 := by
  unfold residue35 residue36 quotient35 modulus
  ring_nf
  all_goals simp
private theorem square36 : residue36 ^ 2 = residue37 + modulus * quotient36 := by
  unfold residue36 residue37 quotient36 modulus
  ring_nf
  all_goals simp
private theorem square37 : residue37 ^ 2 = residue38 + modulus * quotient37 := by
  unfold residue37 residue38 quotient37 modulus
  ring_nf
  all_goals simp
private theorem square38 : residue38 ^ 2 = residue39 + modulus * quotient38 := by
  unfold residue38 residue39 quotient38 modulus
  ring_nf
  all_goals simp
private theorem square39 : residue39 ^ 2 = residue40 + modulus * quotient39 := by
  unfold residue39 residue40 quotient39 modulus
  ring_nf
  all_goals simp
private theorem square40 : residue40 ^ 2 = residue41 + modulus * quotient40 := by
  unfold residue40 residue41 quotient40 modulus
  ring_nf
  all_goals simp
private theorem square41 : residue41 ^ 2 = residue42 + modulus * quotient41 := by
  unfold residue41 residue42 quotient41 modulus
  ring_nf
  all_goals simp
private theorem square42 : residue42 ^ 2 = residue43 + modulus * quotient42 := by
  unfold residue42 residue43 quotient42 modulus
  ring_nf
  all_goals simp
private theorem square43 : residue43 ^ 2 = residue44 + modulus * quotient43 := by
  unfold residue43 residue44 quotient43 modulus
  ring_nf
  all_goals simp
private theorem square44 : residue44 ^ 2 = residue45 + modulus * quotient44 := by
  unfold residue44 residue45 quotient44 modulus
  ring_nf
  all_goals simp
private theorem square45 : residue45 ^ 2 = residue46 + modulus * quotient45 := by
  unfold residue45 residue46 quotient45 modulus
  ring_nf
  all_goals simp
private theorem square46 : residue46 ^ 2 = residue47 + modulus * quotient46 := by
  unfold residue46 residue47 quotient46 modulus
  ring_nf
  all_goals simp
private theorem square47 : residue47 ^ 2 = residue48 + modulus * quotient47 := by
  unfold residue47 residue48 quotient47 modulus
  ring_nf
  all_goals simp
private theorem square48 : residue48 ^ 2 = residue49 + modulus * quotient48 := by
  unfold residue48 residue49 quotient48 modulus
  ring_nf
  all_goals simp
private theorem square49 : residue49 ^ 2 = residue50 + modulus * quotient49 := by
  unfold residue49 residue50 quotient49 modulus
  ring_nf
  all_goals simp
private theorem square50 : residue50 ^ 2 = residue51 + modulus * quotient50 := by
  unfold residue50 residue51 quotient50 modulus
  ring_nf
  all_goals simp
private theorem square51 : residue51 ^ 2 = residue52 + modulus * quotient51 := by
  unfold residue51 residue52 quotient51 modulus
  ring_nf
  all_goals simp
private theorem square52 : residue52 ^ 2 = residue53 + modulus * quotient52 := by
  unfold residue52 residue53 quotient52 modulus
  ring_nf
  all_goals simp
private theorem square53 : residue53 ^ 2 = residue54 + modulus * quotient53 := by
  unfold residue53 residue54 quotient53 modulus
  ring_nf
  all_goals simp
private theorem square54 : residue54 ^ 2 = residue55 + modulus * quotient54 := by
  unfold residue54 residue55 quotient54 modulus
  ring_nf
  all_goals simp
private theorem square55 : residue55 ^ 2 = residue56 + modulus * quotient55 := by
  unfold residue55 residue56 quotient55 modulus
  ring_nf
  all_goals simp
private theorem square56 : residue56 ^ 2 = residue57 + modulus * quotient56 := by
  unfold residue56 residue57 quotient56 modulus
  ring_nf
  all_goals simp
private theorem square57 : residue57 ^ 2 = residue58 + modulus * quotient57 := by
  unfold residue57 residue58 quotient57 modulus
  ring_nf
  all_goals simp
private theorem square58 : residue58 ^ 2 = residue59 + modulus * quotient58 := by
  unfold residue58 residue59 quotient58 modulus
  ring_nf
  all_goals simp
private theorem square59 : residue59 ^ 2 = residue60 + modulus * quotient59 := by
  unfold residue59 residue60 quotient59 modulus
  ring_nf
  all_goals simp
private theorem square60 : residue60 ^ 2 = residue61 + modulus * quotient60 := by
  unfold residue60 residue61 quotient60 modulus
  ring_nf
  all_goals simp
private theorem square61 : residue61 ^ 2 = residue62 + modulus * quotient61 := by
  unfold residue61 residue62 quotient61 modulus
  ring_nf
  all_goals simp
private theorem square62 : residue62 ^ 2 = residue63 + modulus * quotient62 := by
  unfold residue62 residue63 quotient62 modulus
  ring_nf
  all_goals simp
private theorem square63 : residue63 ^ 2 = residue64 + modulus * quotient63 := by
  unfold residue63 residue64 quotient63 modulus
  ring_nf
  all_goals simp
private theorem square64 : residue64 ^ 2 = residue65 + modulus * quotient64 := by
  unfold residue64 residue65 quotient64 modulus
  ring_nf
  all_goals simp
private theorem square65 : residue65 ^ 2 = residue66 + modulus * quotient65 := by
  unfold residue65 residue66 quotient65 modulus
  ring_nf
  all_goals simp
private theorem square66 : residue66 ^ 2 = residue67 + modulus * quotient66 := by
  unfold residue66 residue67 quotient66 modulus
  ring_nf
  all_goals simp
private theorem square67 : residue67 ^ 2 = residue68 + modulus * quotient67 := by
  unfold residue67 residue68 quotient67 modulus
  ring_nf
  all_goals simp
private theorem square68 : residue68 ^ 2 = residue69 + modulus * quotient68 := by
  unfold residue68 residue69 quotient68 modulus
  ring_nf
  all_goals simp
private theorem square69 : residue69 ^ 2 = residue70 + modulus * quotient69 := by
  unfold residue69 residue70 quotient69 modulus
  ring_nf
  all_goals simp
private theorem square70 : residue70 ^ 2 = residue71 + modulus * quotient70 := by
  unfold residue70 residue71 quotient70 modulus
  ring_nf
  all_goals simp
private theorem square71 : residue71 ^ 2 = residue72 + modulus * quotient71 := by
  unfold residue71 residue72 quotient71 modulus
  ring_nf
  all_goals simp
private theorem square72 : residue72 ^ 2 = residue73 + modulus * quotient72 := by
  unfold residue72 residue73 quotient72 modulus
  ring_nf
  all_goals simp
private theorem square73 : residue73 ^ 2 = residue74 + modulus * quotient73 := by
  unfold residue73 residue74 quotient73 modulus
  ring_nf
  all_goals simp
private theorem square74 : residue74 ^ 2 = residue75 + modulus * quotient74 := by
  unfold residue74 residue75 quotient74 modulus
  ring_nf
  all_goals simp
private theorem square75 : residue75 ^ 2 = residue76 + modulus * quotient75 := by
  unfold residue75 residue76 quotient75 modulus
  ring_nf
  all_goals simp
private theorem square76 : residue76 ^ 2 = residue77 + modulus * quotient76 := by
  unfold residue76 residue77 quotient76 modulus
  ring_nf
  all_goals simp
private theorem square77 : residue77 ^ 2 = residue78 + modulus * quotient77 := by
  unfold residue77 residue78 quotient77 modulus
  ring_nf
  all_goals simp
private theorem square78 : residue78 ^ 2 = residue79 + modulus * quotient78 := by
  unfold residue78 residue79 quotient78 modulus
  ring_nf
  all_goals simp
private theorem square79 : residue79 ^ 2 = residue80 + modulus * quotient79 := by
  unfold residue79 residue80 quotient79 modulus
  ring_nf
  all_goals simp
private theorem square80 : residue80 ^ 2 = residue81 + modulus * quotient80 := by
  unfold residue80 residue81 quotient80 modulus
  ring_nf
  all_goals simp
private theorem square81 : residue81 ^ 2 = residue82 + modulus * quotient81 := by
  unfold residue81 residue82 quotient81 modulus
  ring_nf
  all_goals simp
private theorem square82 : residue82 ^ 2 = residue83 + modulus * quotient82 := by
  unfold residue82 residue83 quotient82 modulus
  ring_nf
  all_goals simp
private theorem square83 : residue83 ^ 2 = residue84 + modulus * quotient83 := by
  unfold residue83 residue84 quotient83 modulus
  ring_nf
  all_goals simp
private theorem square84 : residue84 ^ 2 = residue85 + modulus * quotient84 := by
  unfold residue84 residue85 quotient84 modulus
  ring_nf
  all_goals simp
private theorem square85 : residue85 ^ 2 = residue86 + modulus * quotient85 := by
  unfold residue85 residue86 quotient85 modulus
  ring_nf
  all_goals simp
private theorem square86 : residue86 ^ 2 = residue87 + modulus * quotient86 := by
  unfold residue86 residue87 quotient86 modulus
  ring_nf
  all_goals simp
private theorem square87 : residue87 ^ 2 = residue88 + modulus * quotient87 := by
  unfold residue87 residue88 quotient87 modulus
  ring_nf
  all_goals simp
private theorem square88 : residue88 ^ 2 = residue89 + modulus * quotient88 := by
  unfold residue88 residue89 quotient88 modulus
  ring_nf
  all_goals simp
private theorem square89 : residue89 ^ 2 = residue90 + modulus * quotient89 := by
  unfold residue89 residue90 quotient89 modulus
  ring_nf
  all_goals simp
private theorem square90 : residue90 ^ 2 = residue91 + modulus * quotient90 := by
  unfold residue90 residue91 quotient90 modulus
  ring_nf
  all_goals simp
private theorem square91 : residue91 ^ 2 = residue92 + modulus * quotient91 := by
  unfold residue91 residue92 quotient91 modulus
  ring_nf
  all_goals simp
private theorem square92 : residue92 ^ 2 = residue93 + modulus * quotient92 := by
  unfold residue92 residue93 quotient92 modulus
  ring_nf
  all_goals simp
private theorem square93 : residue93 ^ 2 = residue94 + modulus * quotient93 := by
  unfold residue93 residue94 quotient93 modulus
  ring_nf
  all_goals simp
private theorem square94 : residue94 ^ 2 = residue95 + modulus * quotient94 := by
  unfold residue94 residue95 quotient94 modulus
  ring_nf
  all_goals simp
private theorem square95 : residue95 ^ 2 = residue96 + modulus * quotient95 := by
  unfold residue95 residue96 quotient95 modulus
  ring_nf
  all_goals simp
private theorem square96 : residue96 ^ 2 = residue97 + modulus * quotient96 := by
  unfold residue96 residue97 quotient96 modulus
  ring_nf
  all_goals simp
private theorem square97 : residue97 ^ 2 = residue98 + modulus * quotient97 := by
  unfold residue97 residue98 quotient97 modulus
  ring_nf
  all_goals simp
private theorem square98 : residue98 ^ 2 = residue99 + modulus * quotient98 := by
  unfold residue98 residue99 quotient98 modulus
  ring_nf
  all_goals simp
private theorem square99 : residue99 ^ 2 = residue100 + modulus * quotient99 := by
  unfold residue99 residue100 quotient99 modulus
  ring_nf
  all_goals simp
private theorem square100 : residue100 ^ 2 = residue101 + modulus * quotient100 := by
  unfold residue100 residue101 quotient100 modulus
  ring_nf
  all_goals simp
private theorem square101 : residue101 ^ 2 = residue102 + modulus * quotient101 := by
  unfold residue101 residue102 quotient101 modulus
  ring_nf
  all_goals simp
private theorem square102 : residue102 ^ 2 = residue103 + modulus * quotient102 := by
  unfold residue102 residue103 quotient102 modulus
  ring_nf
  all_goals simp
private theorem square103 : residue103 ^ 2 = residue104 + modulus * quotient103 := by
  unfold residue103 residue104 quotient103 modulus
  ring_nf
  all_goals simp
private theorem square104 : residue104 ^ 2 = residue105 + modulus * quotient104 := by
  unfold residue104 residue105 quotient104 modulus
  ring_nf
  all_goals simp
private theorem square105 : residue105 ^ 2 = residue106 + modulus * quotient105 := by
  unfold residue105 residue106 quotient105 modulus
  ring_nf
  all_goals simp
private theorem square106 : residue106 ^ 2 = residue107 + modulus * quotient106 := by
  unfold residue106 residue107 quotient106 modulus
  ring_nf
  all_goals simp
private theorem square107 : residue107 ^ 2 = residue108 + modulus * quotient107 := by
  unfold residue107 residue108 quotient107 modulus
  ring_nf
  all_goals simp
private theorem square108 : residue108 ^ 2 = residue109 + modulus * quotient108 := by
  unfold residue108 residue109 quotient108 modulus
  ring_nf
  all_goals simp
private theorem square109 : residue109 ^ 2 = residue110 + modulus * quotient109 := by
  unfold residue109 residue110 quotient109 modulus
  ring_nf
  all_goals simp
private theorem square110 : residue110 ^ 2 = residue111 + modulus * quotient110 := by
  unfold residue110 residue111 quotient110 modulus
  ring_nf
  all_goals simp
private theorem square111 : residue111 ^ 2 = residue112 + modulus * quotient111 := by
  unfold residue111 residue112 quotient111 modulus
  ring_nf
  all_goals simp
private theorem square112 : residue112 ^ 2 = residue113 + modulus * quotient112 := by
  unfold residue112 residue113 quotient112 modulus
  ring_nf
  all_goals simp
private theorem square113 : residue113 ^ 2 = residue114 + modulus * quotient113 := by
  unfold residue113 residue114 quotient113 modulus
  ring_nf
  all_goals simp
private theorem square114 : residue114 ^ 2 = residue115 + modulus * quotient114 := by
  unfold residue114 residue115 quotient114 modulus
  ring_nf
  all_goals simp
private theorem square115 : residue115 ^ 2 = residue116 + modulus * quotient115 := by
  unfold residue115 residue116 quotient115 modulus
  ring_nf
  all_goals simp
private theorem square116 : residue116 ^ 2 = residue117 + modulus * quotient116 := by
  unfold residue116 residue117 quotient116 modulus
  ring_nf
  all_goals simp
private theorem square117 : residue117 ^ 2 = residue118 + modulus * quotient117 := by
  unfold residue117 residue118 quotient117 modulus
  ring_nf
  all_goals simp
private theorem square118 : residue118 ^ 2 = residue119 + modulus * quotient118 := by
  unfold residue118 residue119 quotient118 modulus
  ring_nf
  all_goals simp
private theorem square119 : residue119 ^ 2 = residue120 + modulus * quotient119 := by
  unfold residue119 residue120 quotient119 modulus
  ring_nf
  all_goals simp
private theorem square120 : residue120 ^ 2 = residue121 + modulus * quotient120 := by
  unfold residue120 residue121 quotient120 modulus
  ring_nf
  all_goals simp
private theorem square121 : residue121 ^ 2 = residue122 + modulus * quotient121 := by
  unfold residue121 residue122 quotient121 modulus
  ring_nf
  all_goals simp
private theorem square122 : residue122 ^ 2 = residue123 + modulus * quotient122 := by
  unfold residue122 residue123 quotient122 modulus
  ring_nf
  all_goals simp
private theorem square123 : residue123 ^ 2 = residue124 + modulus * quotient123 := by
  unfold residue123 residue124 quotient123 modulus
  ring_nf
  all_goals simp
private theorem square124 : residue124 ^ 2 = residue125 + modulus * quotient124 := by
  unfold residue124 residue125 quotient124 modulus
  ring_nf
  all_goals simp
private theorem square125 : residue125 ^ 2 = residue126 + modulus * quotient125 := by
  unfold residue125 residue126 quotient125 modulus
  ring_nf
  all_goals simp
private theorem square126 : residue126 ^ 2 = residue127 + modulus * quotient126 := by
  unfold residue126 residue127 quotient126 modulus
  ring_nf
  all_goals simp
private theorem square127 : residue127 ^ 2 = residue128 + modulus * quotient127 := by
  unfold residue127 residue128 quotient127 modulus
  ring_nf
  all_goals simp


private irreducible_def twoPower (i : Nat) : Nat := 2 ^ i

private theorem twoPower_succ (i : Nat) : twoPower (i + 1) = twoPower i * 2 := by
  rw [twoPower_def, twoPower_def, pow_succ]

private theorem frobeniusStep (f r s q : F₂[X]) (i : Nat)
    (h : f ∣ X ^ twoPower i - r) (hsq : r ^ 2 = s + f * q) :
    f ∣ X ^ twoPower (i + 1) - s := by
  obtain ⟨w, hw⟩ := h
  refine ⟨w * (X ^ twoPower i + r) + q, ?_⟩
  have hp : X ^ twoPower (i + 1) = (X ^ twoPower i : F₂[X]) ^ 2 := by
    rw [twoPower_succ, pow_mul]
  calc
    X ^ twoPower (i + 1) - s = (X ^ twoPower i) ^ 2 - s := by rw [hp]
    _ = (X ^ twoPower i - r) * (X ^ twoPower i + r) + f * q := by
      linear_combination hsq
    _ = f * (w * (X ^ twoPower i + r) + q) := by rw [hw]; ring

private theorem frobenius0 : modulus ∣ X ^ twoPower 0 - residue0 := by
  rw [twoPower_def]
  simp [residue0]
private theorem frobenius1 : modulus ∣ X ^ twoPower 1 - residue1 :=
  frobeniusStep modulus residue0 residue1 quotient0 0 frobenius0 square0
private theorem frobenius2 : modulus ∣ X ^ twoPower 2 - residue2 :=
  frobeniusStep modulus residue1 residue2 quotient1 1 frobenius1 square1
private theorem frobenius3 : modulus ∣ X ^ twoPower 3 - residue3 :=
  frobeniusStep modulus residue2 residue3 quotient2 2 frobenius2 square2
private theorem frobenius4 : modulus ∣ X ^ twoPower 4 - residue4 :=
  frobeniusStep modulus residue3 residue4 quotient3 3 frobenius3 square3
private theorem frobenius5 : modulus ∣ X ^ twoPower 5 - residue5 :=
  frobeniusStep modulus residue4 residue5 quotient4 4 frobenius4 square4
private theorem frobenius6 : modulus ∣ X ^ twoPower 6 - residue6 :=
  frobeniusStep modulus residue5 residue6 quotient5 5 frobenius5 square5
private theorem frobenius7 : modulus ∣ X ^ twoPower 7 - residue7 :=
  frobeniusStep modulus residue6 residue7 quotient6 6 frobenius6 square6
private theorem frobenius8 : modulus ∣ X ^ twoPower 8 - residue8 :=
  frobeniusStep modulus residue7 residue8 quotient7 7 frobenius7 square7
private theorem frobenius9 : modulus ∣ X ^ twoPower 9 - residue9 :=
  frobeniusStep modulus residue8 residue9 quotient8 8 frobenius8 square8
private theorem frobenius10 : modulus ∣ X ^ twoPower 10 - residue10 :=
  frobeniusStep modulus residue9 residue10 quotient9 9 frobenius9 square9
private theorem frobenius11 : modulus ∣ X ^ twoPower 11 - residue11 :=
  frobeniusStep modulus residue10 residue11 quotient10 10 frobenius10 square10
private theorem frobenius12 : modulus ∣ X ^ twoPower 12 - residue12 :=
  frobeniusStep modulus residue11 residue12 quotient11 11 frobenius11 square11
private theorem frobenius13 : modulus ∣ X ^ twoPower 13 - residue13 :=
  frobeniusStep modulus residue12 residue13 quotient12 12 frobenius12 square12
private theorem frobenius14 : modulus ∣ X ^ twoPower 14 - residue14 :=
  frobeniusStep modulus residue13 residue14 quotient13 13 frobenius13 square13
private theorem frobenius15 : modulus ∣ X ^ twoPower 15 - residue15 :=
  frobeniusStep modulus residue14 residue15 quotient14 14 frobenius14 square14
private theorem frobenius16 : modulus ∣ X ^ twoPower 16 - residue16 :=
  frobeniusStep modulus residue15 residue16 quotient15 15 frobenius15 square15
private theorem frobenius17 : modulus ∣ X ^ twoPower 17 - residue17 :=
  frobeniusStep modulus residue16 residue17 quotient16 16 frobenius16 square16
private theorem frobenius18 : modulus ∣ X ^ twoPower 18 - residue18 :=
  frobeniusStep modulus residue17 residue18 quotient17 17 frobenius17 square17
private theorem frobenius19 : modulus ∣ X ^ twoPower 19 - residue19 :=
  frobeniusStep modulus residue18 residue19 quotient18 18 frobenius18 square18
private theorem frobenius20 : modulus ∣ X ^ twoPower 20 - residue20 :=
  frobeniusStep modulus residue19 residue20 quotient19 19 frobenius19 square19
private theorem frobenius21 : modulus ∣ X ^ twoPower 21 - residue21 :=
  frobeniusStep modulus residue20 residue21 quotient20 20 frobenius20 square20
private theorem frobenius22 : modulus ∣ X ^ twoPower 22 - residue22 :=
  frobeniusStep modulus residue21 residue22 quotient21 21 frobenius21 square21
private theorem frobenius23 : modulus ∣ X ^ twoPower 23 - residue23 :=
  frobeniusStep modulus residue22 residue23 quotient22 22 frobenius22 square22
private theorem frobenius24 : modulus ∣ X ^ twoPower 24 - residue24 :=
  frobeniusStep modulus residue23 residue24 quotient23 23 frobenius23 square23
private theorem frobenius25 : modulus ∣ X ^ twoPower 25 - residue25 :=
  frobeniusStep modulus residue24 residue25 quotient24 24 frobenius24 square24
private theorem frobenius26 : modulus ∣ X ^ twoPower 26 - residue26 :=
  frobeniusStep modulus residue25 residue26 quotient25 25 frobenius25 square25
private theorem frobenius27 : modulus ∣ X ^ twoPower 27 - residue27 :=
  frobeniusStep modulus residue26 residue27 quotient26 26 frobenius26 square26
private theorem frobenius28 : modulus ∣ X ^ twoPower 28 - residue28 :=
  frobeniusStep modulus residue27 residue28 quotient27 27 frobenius27 square27
private theorem frobenius29 : modulus ∣ X ^ twoPower 29 - residue29 :=
  frobeniusStep modulus residue28 residue29 quotient28 28 frobenius28 square28
private theorem frobenius30 : modulus ∣ X ^ twoPower 30 - residue30 :=
  frobeniusStep modulus residue29 residue30 quotient29 29 frobenius29 square29
private theorem frobenius31 : modulus ∣ X ^ twoPower 31 - residue31 :=
  frobeniusStep modulus residue30 residue31 quotient30 30 frobenius30 square30
private theorem frobenius32 : modulus ∣ X ^ twoPower 32 - residue32 :=
  frobeniusStep modulus residue31 residue32 quotient31 31 frobenius31 square31
private theorem frobenius33 : modulus ∣ X ^ twoPower 33 - residue33 :=
  frobeniusStep modulus residue32 residue33 quotient32 32 frobenius32 square32
private theorem frobenius34 : modulus ∣ X ^ twoPower 34 - residue34 :=
  frobeniusStep modulus residue33 residue34 quotient33 33 frobenius33 square33
private theorem frobenius35 : modulus ∣ X ^ twoPower 35 - residue35 :=
  frobeniusStep modulus residue34 residue35 quotient34 34 frobenius34 square34
private theorem frobenius36 : modulus ∣ X ^ twoPower 36 - residue36 :=
  frobeniusStep modulus residue35 residue36 quotient35 35 frobenius35 square35
private theorem frobenius37 : modulus ∣ X ^ twoPower 37 - residue37 :=
  frobeniusStep modulus residue36 residue37 quotient36 36 frobenius36 square36
private theorem frobenius38 : modulus ∣ X ^ twoPower 38 - residue38 :=
  frobeniusStep modulus residue37 residue38 quotient37 37 frobenius37 square37
private theorem frobenius39 : modulus ∣ X ^ twoPower 39 - residue39 :=
  frobeniusStep modulus residue38 residue39 quotient38 38 frobenius38 square38
private theorem frobenius40 : modulus ∣ X ^ twoPower 40 - residue40 :=
  frobeniusStep modulus residue39 residue40 quotient39 39 frobenius39 square39
private theorem frobenius41 : modulus ∣ X ^ twoPower 41 - residue41 :=
  frobeniusStep modulus residue40 residue41 quotient40 40 frobenius40 square40
private theorem frobenius42 : modulus ∣ X ^ twoPower 42 - residue42 :=
  frobeniusStep modulus residue41 residue42 quotient41 41 frobenius41 square41
private theorem frobenius43 : modulus ∣ X ^ twoPower 43 - residue43 :=
  frobeniusStep modulus residue42 residue43 quotient42 42 frobenius42 square42
private theorem frobenius44 : modulus ∣ X ^ twoPower 44 - residue44 :=
  frobeniusStep modulus residue43 residue44 quotient43 43 frobenius43 square43
private theorem frobenius45 : modulus ∣ X ^ twoPower 45 - residue45 :=
  frobeniusStep modulus residue44 residue45 quotient44 44 frobenius44 square44
private theorem frobenius46 : modulus ∣ X ^ twoPower 46 - residue46 :=
  frobeniusStep modulus residue45 residue46 quotient45 45 frobenius45 square45
private theorem frobenius47 : modulus ∣ X ^ twoPower 47 - residue47 :=
  frobeniusStep modulus residue46 residue47 quotient46 46 frobenius46 square46
private theorem frobenius48 : modulus ∣ X ^ twoPower 48 - residue48 :=
  frobeniusStep modulus residue47 residue48 quotient47 47 frobenius47 square47
private theorem frobenius49 : modulus ∣ X ^ twoPower 49 - residue49 :=
  frobeniusStep modulus residue48 residue49 quotient48 48 frobenius48 square48
private theorem frobenius50 : modulus ∣ X ^ twoPower 50 - residue50 :=
  frobeniusStep modulus residue49 residue50 quotient49 49 frobenius49 square49
private theorem frobenius51 : modulus ∣ X ^ twoPower 51 - residue51 :=
  frobeniusStep modulus residue50 residue51 quotient50 50 frobenius50 square50
private theorem frobenius52 : modulus ∣ X ^ twoPower 52 - residue52 :=
  frobeniusStep modulus residue51 residue52 quotient51 51 frobenius51 square51
private theorem frobenius53 : modulus ∣ X ^ twoPower 53 - residue53 :=
  frobeniusStep modulus residue52 residue53 quotient52 52 frobenius52 square52
private theorem frobenius54 : modulus ∣ X ^ twoPower 54 - residue54 :=
  frobeniusStep modulus residue53 residue54 quotient53 53 frobenius53 square53
private theorem frobenius55 : modulus ∣ X ^ twoPower 55 - residue55 :=
  frobeniusStep modulus residue54 residue55 quotient54 54 frobenius54 square54
private theorem frobenius56 : modulus ∣ X ^ twoPower 56 - residue56 :=
  frobeniusStep modulus residue55 residue56 quotient55 55 frobenius55 square55
private theorem frobenius57 : modulus ∣ X ^ twoPower 57 - residue57 :=
  frobeniusStep modulus residue56 residue57 quotient56 56 frobenius56 square56
private theorem frobenius58 : modulus ∣ X ^ twoPower 58 - residue58 :=
  frobeniusStep modulus residue57 residue58 quotient57 57 frobenius57 square57
private theorem frobenius59 : modulus ∣ X ^ twoPower 59 - residue59 :=
  frobeniusStep modulus residue58 residue59 quotient58 58 frobenius58 square58
private theorem frobenius60 : modulus ∣ X ^ twoPower 60 - residue60 :=
  frobeniusStep modulus residue59 residue60 quotient59 59 frobenius59 square59
private theorem frobenius61 : modulus ∣ X ^ twoPower 61 - residue61 :=
  frobeniusStep modulus residue60 residue61 quotient60 60 frobenius60 square60
private theorem frobenius62 : modulus ∣ X ^ twoPower 62 - residue62 :=
  frobeniusStep modulus residue61 residue62 quotient61 61 frobenius61 square61
private theorem frobenius63 : modulus ∣ X ^ twoPower 63 - residue63 :=
  frobeniusStep modulus residue62 residue63 quotient62 62 frobenius62 square62
private theorem frobenius64 : modulus ∣ X ^ twoPower 64 - residue64 :=
  frobeniusStep modulus residue63 residue64 quotient63 63 frobenius63 square63
private theorem frobenius65 : modulus ∣ X ^ twoPower 65 - residue65 :=
  frobeniusStep modulus residue64 residue65 quotient64 64 frobenius64 square64
private theorem frobenius66 : modulus ∣ X ^ twoPower 66 - residue66 :=
  frobeniusStep modulus residue65 residue66 quotient65 65 frobenius65 square65
private theorem frobenius67 : modulus ∣ X ^ twoPower 67 - residue67 :=
  frobeniusStep modulus residue66 residue67 quotient66 66 frobenius66 square66
private theorem frobenius68 : modulus ∣ X ^ twoPower 68 - residue68 :=
  frobeniusStep modulus residue67 residue68 quotient67 67 frobenius67 square67
private theorem frobenius69 : modulus ∣ X ^ twoPower 69 - residue69 :=
  frobeniusStep modulus residue68 residue69 quotient68 68 frobenius68 square68
private theorem frobenius70 : modulus ∣ X ^ twoPower 70 - residue70 :=
  frobeniusStep modulus residue69 residue70 quotient69 69 frobenius69 square69
private theorem frobenius71 : modulus ∣ X ^ twoPower 71 - residue71 :=
  frobeniusStep modulus residue70 residue71 quotient70 70 frobenius70 square70
private theorem frobenius72 : modulus ∣ X ^ twoPower 72 - residue72 :=
  frobeniusStep modulus residue71 residue72 quotient71 71 frobenius71 square71
private theorem frobenius73 : modulus ∣ X ^ twoPower 73 - residue73 :=
  frobeniusStep modulus residue72 residue73 quotient72 72 frobenius72 square72
private theorem frobenius74 : modulus ∣ X ^ twoPower 74 - residue74 :=
  frobeniusStep modulus residue73 residue74 quotient73 73 frobenius73 square73
private theorem frobenius75 : modulus ∣ X ^ twoPower 75 - residue75 :=
  frobeniusStep modulus residue74 residue75 quotient74 74 frobenius74 square74
private theorem frobenius76 : modulus ∣ X ^ twoPower 76 - residue76 :=
  frobeniusStep modulus residue75 residue76 quotient75 75 frobenius75 square75
private theorem frobenius77 : modulus ∣ X ^ twoPower 77 - residue77 :=
  frobeniusStep modulus residue76 residue77 quotient76 76 frobenius76 square76
private theorem frobenius78 : modulus ∣ X ^ twoPower 78 - residue78 :=
  frobeniusStep modulus residue77 residue78 quotient77 77 frobenius77 square77
private theorem frobenius79 : modulus ∣ X ^ twoPower 79 - residue79 :=
  frobeniusStep modulus residue78 residue79 quotient78 78 frobenius78 square78
private theorem frobenius80 : modulus ∣ X ^ twoPower 80 - residue80 :=
  frobeniusStep modulus residue79 residue80 quotient79 79 frobenius79 square79
private theorem frobenius81 : modulus ∣ X ^ twoPower 81 - residue81 :=
  frobeniusStep modulus residue80 residue81 quotient80 80 frobenius80 square80
private theorem frobenius82 : modulus ∣ X ^ twoPower 82 - residue82 :=
  frobeniusStep modulus residue81 residue82 quotient81 81 frobenius81 square81
private theorem frobenius83 : modulus ∣ X ^ twoPower 83 - residue83 :=
  frobeniusStep modulus residue82 residue83 quotient82 82 frobenius82 square82
private theorem frobenius84 : modulus ∣ X ^ twoPower 84 - residue84 :=
  frobeniusStep modulus residue83 residue84 quotient83 83 frobenius83 square83
private theorem frobenius85 : modulus ∣ X ^ twoPower 85 - residue85 :=
  frobeniusStep modulus residue84 residue85 quotient84 84 frobenius84 square84
private theorem frobenius86 : modulus ∣ X ^ twoPower 86 - residue86 :=
  frobeniusStep modulus residue85 residue86 quotient85 85 frobenius85 square85
private theorem frobenius87 : modulus ∣ X ^ twoPower 87 - residue87 :=
  frobeniusStep modulus residue86 residue87 quotient86 86 frobenius86 square86
private theorem frobenius88 : modulus ∣ X ^ twoPower 88 - residue88 :=
  frobeniusStep modulus residue87 residue88 quotient87 87 frobenius87 square87
private theorem frobenius89 : modulus ∣ X ^ twoPower 89 - residue89 :=
  frobeniusStep modulus residue88 residue89 quotient88 88 frobenius88 square88
private theorem frobenius90 : modulus ∣ X ^ twoPower 90 - residue90 :=
  frobeniusStep modulus residue89 residue90 quotient89 89 frobenius89 square89
private theorem frobenius91 : modulus ∣ X ^ twoPower 91 - residue91 :=
  frobeniusStep modulus residue90 residue91 quotient90 90 frobenius90 square90
private theorem frobenius92 : modulus ∣ X ^ twoPower 92 - residue92 :=
  frobeniusStep modulus residue91 residue92 quotient91 91 frobenius91 square91
private theorem frobenius93 : modulus ∣ X ^ twoPower 93 - residue93 :=
  frobeniusStep modulus residue92 residue93 quotient92 92 frobenius92 square92
private theorem frobenius94 : modulus ∣ X ^ twoPower 94 - residue94 :=
  frobeniusStep modulus residue93 residue94 quotient93 93 frobenius93 square93
private theorem frobenius95 : modulus ∣ X ^ twoPower 95 - residue95 :=
  frobeniusStep modulus residue94 residue95 quotient94 94 frobenius94 square94
private theorem frobenius96 : modulus ∣ X ^ twoPower 96 - residue96 :=
  frobeniusStep modulus residue95 residue96 quotient95 95 frobenius95 square95
private theorem frobenius97 : modulus ∣ X ^ twoPower 97 - residue97 :=
  frobeniusStep modulus residue96 residue97 quotient96 96 frobenius96 square96
private theorem frobenius98 : modulus ∣ X ^ twoPower 98 - residue98 :=
  frobeniusStep modulus residue97 residue98 quotient97 97 frobenius97 square97
private theorem frobenius99 : modulus ∣ X ^ twoPower 99 - residue99 :=
  frobeniusStep modulus residue98 residue99 quotient98 98 frobenius98 square98
private theorem frobenius100 : modulus ∣ X ^ twoPower 100 - residue100 :=
  frobeniusStep modulus residue99 residue100 quotient99 99 frobenius99 square99
private theorem frobenius101 : modulus ∣ X ^ twoPower 101 - residue101 :=
  frobeniusStep modulus residue100 residue101 quotient100 100 frobenius100 square100
private theorem frobenius102 : modulus ∣ X ^ twoPower 102 - residue102 :=
  frobeniusStep modulus residue101 residue102 quotient101 101 frobenius101 square101
private theorem frobenius103 : modulus ∣ X ^ twoPower 103 - residue103 :=
  frobeniusStep modulus residue102 residue103 quotient102 102 frobenius102 square102
private theorem frobenius104 : modulus ∣ X ^ twoPower 104 - residue104 :=
  frobeniusStep modulus residue103 residue104 quotient103 103 frobenius103 square103
private theorem frobenius105 : modulus ∣ X ^ twoPower 105 - residue105 :=
  frobeniusStep modulus residue104 residue105 quotient104 104 frobenius104 square104
private theorem frobenius106 : modulus ∣ X ^ twoPower 106 - residue106 :=
  frobeniusStep modulus residue105 residue106 quotient105 105 frobenius105 square105
private theorem frobenius107 : modulus ∣ X ^ twoPower 107 - residue107 :=
  frobeniusStep modulus residue106 residue107 quotient106 106 frobenius106 square106
private theorem frobenius108 : modulus ∣ X ^ twoPower 108 - residue108 :=
  frobeniusStep modulus residue107 residue108 quotient107 107 frobenius107 square107
private theorem frobenius109 : modulus ∣ X ^ twoPower 109 - residue109 :=
  frobeniusStep modulus residue108 residue109 quotient108 108 frobenius108 square108
private theorem frobenius110 : modulus ∣ X ^ twoPower 110 - residue110 :=
  frobeniusStep modulus residue109 residue110 quotient109 109 frobenius109 square109
private theorem frobenius111 : modulus ∣ X ^ twoPower 111 - residue111 :=
  frobeniusStep modulus residue110 residue111 quotient110 110 frobenius110 square110
private theorem frobenius112 : modulus ∣ X ^ twoPower 112 - residue112 :=
  frobeniusStep modulus residue111 residue112 quotient111 111 frobenius111 square111
private theorem frobenius113 : modulus ∣ X ^ twoPower 113 - residue113 :=
  frobeniusStep modulus residue112 residue113 quotient112 112 frobenius112 square112
private theorem frobenius114 : modulus ∣ X ^ twoPower 114 - residue114 :=
  frobeniusStep modulus residue113 residue114 quotient113 113 frobenius113 square113
private theorem frobenius115 : modulus ∣ X ^ twoPower 115 - residue115 :=
  frobeniusStep modulus residue114 residue115 quotient114 114 frobenius114 square114
private theorem frobenius116 : modulus ∣ X ^ twoPower 116 - residue116 :=
  frobeniusStep modulus residue115 residue116 quotient115 115 frobenius115 square115
private theorem frobenius117 : modulus ∣ X ^ twoPower 117 - residue117 :=
  frobeniusStep modulus residue116 residue117 quotient116 116 frobenius116 square116
private theorem frobenius118 : modulus ∣ X ^ twoPower 118 - residue118 :=
  frobeniusStep modulus residue117 residue118 quotient117 117 frobenius117 square117
private theorem frobenius119 : modulus ∣ X ^ twoPower 119 - residue119 :=
  frobeniusStep modulus residue118 residue119 quotient118 118 frobenius118 square118
private theorem frobenius120 : modulus ∣ X ^ twoPower 120 - residue120 :=
  frobeniusStep modulus residue119 residue120 quotient119 119 frobenius119 square119
private theorem frobenius121 : modulus ∣ X ^ twoPower 121 - residue121 :=
  frobeniusStep modulus residue120 residue121 quotient120 120 frobenius120 square120
private theorem frobenius122 : modulus ∣ X ^ twoPower 122 - residue122 :=
  frobeniusStep modulus residue121 residue122 quotient121 121 frobenius121 square121
private theorem frobenius123 : modulus ∣ X ^ twoPower 123 - residue123 :=
  frobeniusStep modulus residue122 residue123 quotient122 122 frobenius122 square122
private theorem frobenius124 : modulus ∣ X ^ twoPower 124 - residue124 :=
  frobeniusStep modulus residue123 residue124 quotient123 123 frobenius123 square123
private theorem frobenius125 : modulus ∣ X ^ twoPower 125 - residue125 :=
  frobeniusStep modulus residue124 residue125 quotient124 124 frobenius124 square124
private theorem frobenius126 : modulus ∣ X ^ twoPower 126 - residue126 :=
  frobeniusStep modulus residue125 residue126 quotient125 125 frobenius125 square125
private theorem frobenius127 : modulus ∣ X ^ twoPower 127 - residue127 :=
  frobeniusStep modulus residue126 residue127 quotient126 126 frobenius126 square126
private theorem frobenius128 : modulus ∣ X ^ twoPower 128 - residue128 :=
  frobeniusStep modulus residue127 residue128 quotient127 127 frobenius127 square127

private theorem frobeniusDivisibility :
    (modulus ∣ X ^ twoPower 64 - residue64) ∧
      (modulus ∣ X ^ twoPower 128 - residue128) :=
  ⟨frobenius64, frobenius128⟩

private noncomputable def bezoutLeft : F₂[X] := 1 + X + X ^ 3 + X ^ 4 + X ^ 5 + X ^ 6 + X ^ 10 + X ^ 12 + X ^ 14 + X ^ 16 + X ^ 17 + X ^ 21 + X ^ 22 + X ^ 23 + X ^ 24 + X ^ 25 + X ^ 28 + X ^ 29 + X ^ 30 + X ^ 32 + X ^ 33 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 40 + X ^ 42 + X ^ 46 + X ^ 47 + X ^ 49 + X ^ 50 + X ^ 51 + X ^ 54 + X ^ 56 + X ^ 61 + X ^ 65 + X ^ 66 + X ^ 68 + X ^ 69 + X ^ 75 + X ^ 76 + X ^ 78 + X ^ 83 + X ^ 85 + X ^ 86 + X ^ 88 + X ^ 89 + X ^ 91 + X ^ 92 + X ^ 93 + X ^ 94 + X ^ 95 + X ^ 96 + X ^ 98 + X ^ 101 + X ^ 103 + X ^ 104 + X ^ 107 + X ^ 108 + X ^ 112 + X ^ 113 + X ^ 116 + X ^ 117 + X ^ 118 + X ^ 119 + X ^ 122 + X ^ 124 + X ^ 125
private noncomputable def bezoutRight : F₂[X] := X ^ 3 + X ^ 4 + X ^ 12 + X ^ 13 + X ^ 15 + X ^ 16 + X ^ 18 + X ^ 20 + X ^ 23 + X ^ 26 + X ^ 27 + X ^ 29 + X ^ 32 + X ^ 34 + X ^ 35 + X ^ 36 + X ^ 37 + X ^ 38 + X ^ 39 + X ^ 41 + X ^ 43 + X ^ 46 + X ^ 47 + X ^ 51 + X ^ 53 + X ^ 55 + X ^ 56 + X ^ 62 + X ^ 65 + X ^ 67 + X ^ 71 + X ^ 72 + X ^ 73 + X ^ 79 + X ^ 82 + X ^ 83 + X ^ 84 + X ^ 85 + X ^ 87 + X ^ 88 + X ^ 90 + X ^ 91 + X ^ 92 + X ^ 94 + X ^ 95 + X ^ 99 + X ^ 101 + X ^ 104 + X ^ 107 + X ^ 108 + X ^ 109 + X ^ 110 + X ^ 111 + X ^ 113 + X ^ 115 + X ^ 118 + X ^ 121 + X ^ 122 + X ^ 123 + X ^ 124 + X ^ 127

private theorem residue64Coprime : IsCoprime modulus (residue64 - X) := by
  refine ⟨bezoutLeft, bezoutRight, ?_⟩
  unfold bezoutLeft bezoutRight residue64 modulus
  ring_nf
  simp

private theorem coprimeLift (f r : F₂[X]) (i : Nat)
    (hc : IsCoprime f (r - X)) (h : f ∣ X ^ twoPower i - r) :
    IsCoprime f (X ^ twoPower i - X) := by
  obtain ⟨w, hw⟩ := h
  have heq : X ^ twoPower i - X = (r - X) + w * f := by
    linear_combination hw
  rw [heq]
  exact hc.add_mul_right_right w

theorem modulus_splits : modulus ∣ X ^ (Nat.card F₂) ^ 128 - X := by
  simpa [F₂, residue128, twoPower_def] using frobeniusDivisibility.2

theorem modulus_coprime : IsCoprime modulus (X ^ (Nat.card F₂) ^ 64 - X) := by
  simpa [F₂, twoPower_def] using coprimeLift modulus residue64 64 residue64Coprime
    frobeniusDivisibility.1

theorem modulus_irreducible : Irreducible modulus :=
  irreducible_of_rabin_128 modulus_monic modulus_natDegree modulus_splits modulus_coprime

end TzapLean.GF128Proof

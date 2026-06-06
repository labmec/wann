SetFactory("OpenCASCADE"); // Geometry kernel

// Shared parameters used by both geometry generators.
Geometry.MatchMeshTolerance = 1e-08;
Geometry.Tolerance = 1e-10;
Mesh.ToleranceReferenceElement = 1e-10;
General.Terminal = 1;
Geometry.SnapPoints = 1;

DefineConstant[
	Lw = {200.0, Name "Lw"},
    Hw = {10.0, Name "Hw"},
	Rw = {5, Name "Rw"},
	ecc = {0.0, Name "ecc"}, // Eccentricity of the wellbore (-0.5 <= ecc <= 0.5)

	Hr = {20.0, Name "Hr"},
	Wr = {200.0, Name "Wr"},
	Lr = {400.0, Name "Lr"},

	h_div = {3, Name "h_div"}, // Division in z
	axial_div = {31, Name "axial_div"}, // Divisions along the well length

	radial_div = {6, Name "radial_div"}, // Radial division of the near-well region
	p_res = {1.6, Name "p_res"}, // Progression coefficient of the near-well region mesh
	p_well = {0.3, Name "p_well"}, // Progression coefficient of the wellbore mesh

	size_min_res = {15.0, Name "size_min_res"},
	size_max_res = {300.0, Name "size_max_res"},
	dist_min_res = {80.0, Name "dist_min_res"},
	dist_max_res = {200.0, Name "dist_max_res"}
];

// Derived parameters
axial_div = Floor(h_div*Lw/(3*Hr));
size_min_res = Lw/axial_div;
size_max_res = Lr/10;
dist_min_res = Hr;
dist_max_res = Wr/3;

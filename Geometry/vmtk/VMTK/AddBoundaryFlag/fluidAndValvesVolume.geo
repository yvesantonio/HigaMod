Mesh.RemeshAlgorithm = 1; // automatic
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.5;
Mesh.Algorithm = 2; // automatic
Mesh.Algorithm3D = 1; // automatic
Mesh.Optimize = 1;
Mesh.RemeshParametrization = 7;
Mesh.Smoothing = 10;

Merge "fluidAndValves.msh";
Geometry.Tolerance=1e-6;
Coherence Mesh;

Surface Loop(16)={13,-14,-15,3,2};
Volume(17)={16};

Physical Volume(18)={17};

Mesh.RemeshAlgorithm = 1; // automatic
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.5;
Mesh.Algorithm = 2; // automatic
Mesh.Algorithm3D = 4; // automatic
Mesh.Optimize = 1;
Mesh.RemeshParametrization = 7;
Mesh.Smoothing = 10;

Merge "fluid.msh";
Merge "surface2.stl";
Merge "surface3.stl";
CreateTopology;

Compound Line(8)={6};
Compound Line(9)={7};

Line Loop(10)={8};
Line Loop(11)={9};
Line Loop(12)={5};

Plane Surface(13)={12,10,11};
Plane Surface(14)={10};
Plane Surface(15)={11};

Physical Surface(16)={13};
Physical Surface(17)={14};
Physical Surface(18)={15};
Physical Surface(19)={3,2};
//Surface Loop(16)={13,14,15,3,2};
//Volume(17)={16};

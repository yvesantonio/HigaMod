Mesh.RemeshAlgorithm = 1; // automatic
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.5;
Mesh.Algorithm = 2; // automatic
Mesh.Algorithm3D = 1; // automatic
Mesh.Optimize = 1;
Mesh.RemeshParametrization = 7;
Mesh.Smoothing = 10;

Merge "fluid.msh";
Merge "mitralValve.stl";
Merge "aorticValve.stl";
Merge "valve_extrusion.stl"; 

CreateTopology;

Compound Line(9)={6};
Compound Line(10)={7};
Compound Line(11)={8};

Line Loop(12)={9};
Line Loop(13)={10};
Line Loop(14)={5};

Plane Surface(15)={14,12,13};
Plane Surface(16)={12};
Plane Surface(17)={13};

Compound Surface(18)={7};

Physical Surface(15)={15};
Physical Surface(16)={16};
Physical Surface(17)={17};
Physical Surface(18)={18};
Physical Surface(1)={3,2};


Surface Loop(19)={15,-16,-17,3,2};
Volume(20)={19};

Surface {18} In Volume {20};

Physical Volume(20)={20};
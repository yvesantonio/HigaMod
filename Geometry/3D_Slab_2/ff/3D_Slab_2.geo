//+
Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {2.0, 0.0, 0.0, 1.0};
//+
Point(3) = {2.0, 3.0, 0.0, 1.0};
//+
Point(4) = {0.0, 3.0, 0.0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 4} {
  Surface{1}; Line{4}; Line{1}; Line{2}; Line{3}; 
}
//+
Physical Surface("101") = {26};
//+
Physical Surface("102") = {13};
//+
Physical Surface("103") = {25};
//+
Physical Surface("104") = {21};
//+
Physical Surface("105") = {17};
//+
Physical Surface("106") = {1};
//+
Physical Surface("106") -= {13};
//+
Physical Surface("106") -= {25};
//+
Physical Surface("106") -= {26};
//+
Physical Surface("106") -= {21};
//+
Physical Surface("106") -= {17};
//+
Physical Surface("106") -= {1};
//+
Physical Surface("101") -= {26};
//+
Physical Surface("102") -= {13};
//+
Physical Surface("103") -= {25};
//+
Physical Surface("104") -= {21};
//+
Physical Surface("105") -= {17};
//+
Physical Surface("1", 1) = {26};
//+
Physical Surface("2", 2) = {13};
//+
Physical Surface("3", 3) = {25};
//+
Physical Surface("4", 4) = {21};
//+
Physical Surface("5", 5) = {17};
//+
Physical Surface("6", 6) = {1};
//+
Physical Volume("7", 7) = {1};

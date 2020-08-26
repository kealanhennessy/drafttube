SetFactory("OpenCASCADE");

Point(20001) = {-13.5, -15.0, 0, 1.0};
Point(20002) = {-13.5, 15.0, 0, 1.0};
Point(20003) = {36.0, 15.0, 0, 1.0};
Point(20004) = {36.0, -15.0, 0, 1.0};
Line(11) = {20001, 20002};
Line(12) = {20002, 20003};
Line(13) = {20003, 20004};
Line(14) = {20001, 20004};

Line Loop(3) = {11, 12, 13, 14};

Point(20005) = {-13.5, -15.0, 1, 1.0};
Point(20006) = {-13.5, 15.0, 1, 1.0};
Point(20007) = {36.0, 15.0, 1, 1.0};
Point(20008) = {36.0, -15.0, 1, 1.0};
Line(15) = {20005, 20006};
Line(16) = {20006, 20007};
Line(17) = {20007, 20008};
Line(18) = {20005, 20008};

Line Loop(4) = {15, 16, 17, 18};

Line(19) = {20001, 20005};
Line(20) = {20002, 20006};
Line(21) = {20003, 20007};
Line(22) = {20004, 20008};

Line Loop(5) = {11, 19, 15, 20}; \\ Line Loops from Lines
Line Loop(6) = {19, 14, 22, 18};
Line Loop(7) = {13, 21, 17, 22};
Line Loop(8) = {12, 20, 16, 21};

Plane Surface(1) = {3}; // Line Loops define Plane Surface
Plane Surface(2) = {4};
Plane Surface(3) = {5};
Plane Surface(4) = {6};
Plane Surface(5) = {7};
Plane Surface(6) = {8};

Surface Loop(1) = {1,2,3,4,5,6}; // Plane Surfaces form a Surface Loop
Volume(1) = {1}; // Surface Loop defines Volume

Physical Surface("Check") = {1}; // Plane Surfaces define Physical Surface, required to define Physical Volume from Volume
Physical Surface("Check2") = {2};
Physical Surface("Check3") = {3};
Physical Surface("Check4") = {4};
Physical Surface("Check5") = {5};
Physical Surface("Check6") = {6};

Physical Volume("PlzWork") = {1}; // Volume defines Physical Volume


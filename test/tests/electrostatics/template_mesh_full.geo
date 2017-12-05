SetFactory("OpenCASCADE");

lc = 1.5;
thickness_para = 5;
thickness_ferro = 10;
layers_para = 10;
layers_ferro = 20;

dimx = 20;
dispx = 0;
dimy = 1;
dispy = 29;

//Points for base plane (bottom plane of bottom paraelectric)
Point(1) = {dispx, dispy, 0, lc};
Point(2) = {dispx + dimx, dispy, 0, lc};
Point(3) = {dispx + dimx, dispy + dimy, 0, lc};
Point(4) = {dispx, dispy + dimy, 0, lc};

//Line loop for base plane
Line(101) = {1, 2};
Line(102) = {2, 3};
Line(103) = {3, 4};
Line(104) = {4, 1};
Line Loop(1) = {101, 102, 103, 104};

//Base surface
Plane Surface(1000) = {1};
Physical Surface("para_bottom_surface_bottom") = {1000};

//Extrusion of base plane to get bottom paraelectric
surfaceVector[] = Extrude {0, 0, thickness_para} {
     Surface{1000};
     Layers{layers_para};
     Recombine;
};

//Get top surface of bottom paraelectric (bottom surface of ferroelectric)
Physical Surface("para_bottom_surface_top") = surfaceVector[0];
//Get volume of bottom paralectric
Physical Volume("para_bottom_volume") = surfaceVector[1];

//Get sides of bottom paraelectric
Physical Surface("para_bottom_side_1") = {surfaceVector[2]};
Physical Surface("para_bottom_side_2") = {surfaceVector[3]};
Physical Surface("para_bottom_side_3") = {surfaceVector[4]};
Physical Surface("para_bottom_side_4") = {surfaceVector[5]};

//Extrusion of bottom paraelectric top plane to get ferroelectric
surfaceVector1[] = Extrude {0, 0, thickness_ferro} {
     Surface{1005};
     Layers{layers_ferro};
     Recombine;
};

//Get top surface of ferroelectric (bottom surface of top paraelectric)
Physical Surface("para_top_surface_bottom") = surfaceVector1[0];
//Get volume of ferroelectric
Physical Volume("ferro_volume") = surfaceVector1[1];

//Get sides of bottom ferroelectric
Physical Surface("ferro_side_1") = {surfaceVector1[2]};
Physical Surface("ferro_side_2") = {surfaceVector1[3]};
Physical Surface("ferro_side_3") = {surfaceVector1[4]};
Physical Surface("ferro_side_4") = {surfaceVector1[5]};

//Extrusion of ferroelectric top plane to get top paraelectric
surfaceVector2[] = Extrude {0, 0, thickness_para} {
     Surface{1010};
     Layers{layers_para};
     Recombine;
};

//Get top surface of top paraelectric
Physical Surface("para_top_surface_top") = surfaceVector2[0];
//Get volume of top paraelectric
Physical Volume("para_top_volume") = surfaceVector2[1];

//Get sides of bottom top paraelectric
Physical Surface("para_top_side_1") = {surfaceVector2[2]};
Physical Surface("para_top_side_2") = {surfaceVector2[3]};
Physical Surface("para_top_side_3") = {surfaceVector2[4]};
Physical Surface("para_top_side_4") = {surfaceVector2[5]};



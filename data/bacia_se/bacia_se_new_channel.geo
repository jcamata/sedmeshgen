// Gmsh project created on Fri Sep 21 09:05:19 2012
//
//                               1 +-----------------------------+ 2
//                                 |                             |
//                                 |                             |
//                                 |                             |
//                                 |                             |
//                                 |                             |
//     10       8                6 |                             |
//     +---------------------------+                             |
//  yc |        |                  . (0,0,0)                     | yt
//     +---------------------------+                             |
//     9  xcs  7     xcw         5 |                             |
//                                 |                             |
//                                 |                             |
//                                 |                             |
//                                 |                             |
//                                 |                             |
//                                 +-----------------------------+
//                                 4              xtw             3


// Fator de Escala
SCALE = 0.001;

// Dimensões do domínio
XMIN = 0.0*SCALE;
XMAX = 14900.0*SCALE;
YMIN = 0.0*SCALE;
YMAX = 11900.0*SCALE;

// Fatores de proporcionalidade
Ccomp = 0.75;
Clarg = 0.2;
Tcomp = 2.5;
Tlarg = 2.0;
CompP    = Ccomp/Tcomp;
LargP    = Tlarg/Tcomp;
AltuP    = 0.5/Tcomp;
lockP    = 0.5/Tcomp;

// Centro do canal
XCANAL = 5000.0*SCALE;
YCANAL = YMAX*SCALE;

// Largura do  canal
LCANAL=500.0*SCALE;



// Níveis de Refinamento
rscal  =  8;
Ref0   =  0.02*rscal;
Ref1   =  0.01*rscal;
Ref2   =  0.01*rscal;

angleZ = 3;
angleX = 45;
angleX1 = 90 - angleX;

XT    =  XMAX-XMIN; 
YT    =  YMAX-YMIN;

// Altura do Tanque
H = YT*AltuP;

YCS   =  YT*Ccomp*0.3;

hcs   =  YCS*Sin(angleZ*Pi/180);

a     =  H*Sin(angleZ*Pi/180);
b     =  H*Cos(angleZ*Pi/180);
f     =  0.25*H*Sin(angleZ*Pi/180);
d     =  0.25*H*Cos(angleZ*Pi/180);

c  = Cos(angleX*Pi/180);
s  = Sin(angleX*Pi/180);
 
O7X = XCANAL-LCANAL;
O7Y = YMAX; 

P7X = XCANAL- LCANAL;
P7Y = YMAX  + YCS;

R7X = c*(P7X-O7X) - s*(P7Y-O7Y) + O7X;
R7Y = s*(P7X-O7X) + c*(P7Y-O7Y) + O7Y;

O8X = XCANAL+LCANAL;
O8Y = YMAX; 
P8X = XCANAL+LCANAL;
P8Y = YMAX  + YCS;

R8X = c*(P8X-O8X) - s*(P8Y-O8Y) + O8X;
R8Y = s*(P8X-O8X) + c*(P8Y-O8Y) + O8Y;

//
O15X = XCANAL-LCANAL;
O15Y = YMAX; 

P15X = XCANAL- LCANAL;
P15Y = YMAX  + YCS - a;

R15X = c*(P15X-O15X) - s*(P15Y-O15Y) + O15X;
R15Y = s*(P15X-O15X) + c*(P15Y-O15Y) + O15Y;

O16X = XCANAL+LCANAL;
O16Y = YMAX; 

P16X = XCANAL+LCANAL;
P16Y = YMAX  + YCS - a;

R16X = c*(P16X-O16X) - s*(P16Y-O16Y) + O16X;
R16Y = s*(P16X-O16X) + c*(P16Y-O16Y) + O16Y;


O17X = XCANAL-LCANAL;
O17Y = YMAX; 

P17X = XCANAL-LCANAL;
P17Y = YMAX  + YCS - f;

R17X = c*(P17X-O17X) - s*(P17Y-O17Y) + O17X;
R17Y = s*(P17X-O17X) + c*(P17Y-O17Y) + O17Y;

O18X = XCANAL+LCANAL;
O18Y = YMAX; 

P18X = XCANAL+LCANAL;
P18Y = YMAX  + YCS - f;

R18X = c*(P18X-O18X) - s*(P18Y-O18Y) + O18X;
R18Y = s*(P18X-O18X) + c*(P18Y-O18Y) + O18Y;

// Rotacao ao redor

R7X = c*(R7X-R8X) - s*(R7Y-R8Y)+R8X;
R7Y = s*(R7X-R8X) + c*(R7Y-R8Y)+R8Y;

R17X = c*(R17X-R18X) - s*(R17Y-R18Y)+R18X;
R17Y = s*(R17X-R18X) + c*(R17Y-R18Y)+R18Y;

R15X = c*(R15X-R16X) - s*(R15Y-R16Y)+R16X;
R15Y = s*(R15X-R16X) + c*(R15Y-R16Y)+R16Y;


//bottom
Point(1)  = {    XMIN,     YMAX         ,  0  ,   Ref1};
Point(2)  = {    XMAX,     YMAX         ,  0  ,   Ref1};
Point(3)  = {    XMAX,          YMIN    ,  0  ,   Ref1};
Point(4)  = {    XMIN,          YMIN    ,  0  ,   Ref1};
Point(5)  = {    XCANAL-LCANAL, YMAX    ,  0  ,   Ref2};
Point(6)  = {    XCANAL+LCANAL, YMAX    ,  0  ,   Ref2};
Point(7)  = {    R7X          , R7Y     ,  hcs,   Ref2};
Point(8)  = {    R8X          , R8Y     ,  hcs,   Ref2};


//top
Point(9)   = {    XMIN,     YMAX         ,  H  ,   Ref1};
Point(10)  = {    XMAX,     YMAX         ,  H  ,   Ref1};
Point(11)  = {    XMAX,          YMIN    ,  H  ,   Ref1};
Point(12)  = {    XMIN,          YMIN    ,  H  ,   Ref1};
Point(13)  = {    XCANAL-LCANAL, YMAX    ,  H  ,   Ref2};
Point(14)  = {    XCANAL+LCANAL, YMAX    ,  H  ,   Ref2};
Point(15)  = {    R15X , R15Y , hcs+b,   Ref2};
Point(16)  = {    R16X , R16Y, hcs+b,   Ref2};


//top sediment
Point(17) = {  R17X ,    R17Y ,  hcs+d,  Ref2};
Point(18) = {  R18X ,    R18Y ,  hcs+d,  Ref2};

Line(1) = {8, 7};
Line(2) = {7, 5};
Line(3) = {5, 1};
Line(4) = {1, 4};
Line(5) = {4, 3};
Line(6) = {3, 2};
Line(7) = {2, 6};
Line(8) = {6, 8};
Line(9) = {16, 15};
Line(10) = {15, 13};
Line(11) = {13, 9};
Line(12) = {9, 12};
Line(13) = {12, 11};
Line(14) = {11, 10};
Line(15) = {10, 14};
Line(16) = {14, 13};
Line(17) = {6, 5};
Line(18) = {14, 16};
Line(19) = {16, 18};
Line(20) = {18, 8};
Line(21) = {15, 17};
Line(22) = {17, 7};
Line(23) = {17, 18};
Line(24) = {9, 1};
Line(25) = {13, 5};
Line(26) = {12, 4};
Line(27) = {11, 3};
Line(28) = {10, 2};
Line(29) = {14, 6};
Line Loop(30) = {20, 1, -22, 23};
Plane Surface(31) = {30};
Line Loop(32) = {23, -19, 9, 21};
Plane Surface(33) = {32};
Line Loop(34) = {1, 2, -17, 8};
Plane Surface(35) = {34};
Line Loop(36) = {22, 2, -25, -10, 21};
Plane Surface(37) = {36};
Line Loop(38) = {8, -20, -19, -18, 29};
Plane Surface(39) = {38};
Line Loop(40) = {9, 10, -16, 18};
Plane Surface(41) = {40};
Line Loop(42) = {17, -25, -16, 29};
Plane Surface(43) = {42};
Line Loop(44) = {3, -24, -11, 25};
Plane Surface(45) = {44};
Line Loop(46) = {24, 4, -26, -12};
Plane Surface(47) = {46};
Line Loop(48) = {5, -27, -13, 26};
Plane Surface(49) = {48};
Line Loop(50) = {6, -28, -14, 27};
Plane Surface(51) = {50};
Line Loop(52) = {7, 17, 3, 4, 5, 6};
Plane Surface(53) = {52};
Line Loop(54) = {7, -29, -15, 28};
Plane Surface(55) = {54};
Line Loop(56) = {15, 16, 11, 12, 13, 14};
Plane Surface(57) = {56};
Surface Loop(58) = {57, 55, 53, 45, 47, 49, 51, 43};
Volume(59) = {58};
Surface Loop(60) = {41, 33, 31, 39, 35, 37, 43};
Volume(61) = {60};

Physical Surface("MOVING") = {53};
Physical Surface("NOSLIP") = {35, 37, 39, 45, 55};
Physical Surface("INLET1") = {31};
Physical Surface("INLET2") = {33};
Physical Surface("OUTFLOWL") = {47};
Physical Surface("OUTFLOWR") = {51};
Physical Surface("OUTFLOWB") = {49};
Physical Surface("PNULL") = {41, 57};
Physical Volume("FLUID") = {59, 61};



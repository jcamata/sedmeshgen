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
angleX = 50;
angleX1 = 90 - angleX;

XT    =  XMAX-XMIN; 
YT    =  YMAX-YMIN;

// Altura do Tanque
H = YT*AltuP;

YCS   =  YT*Ccomp*0.5;

hcs   =  YCS*Sin(angleZ*Pi/180);

a     =  H*Sin(angleZ*Pi/180);
b     =  H*Cos(angleZ*Pi/180);
f     =  0.5*H*Sin(angleZ*Pi/180);
d     =  0.5*H*Cos(angleZ*Pi/180);

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


Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};
Line(4) = {2, 3};
Line(5) = {3, 4};
Line(6) = {4, 1};
Line(7) = {5, 7};
Line(8) = {7, 8};
Line(9) = {8, 6};
Line(10) = {9, 13};
Line(11) = {14, 14};
Line(12) = {13, 14};
Line(13) = {14, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};
Line(17) = {13, 15};
Line(18) = {15, 16};
Line(19) = {16, 14};
Line(20) = {9, 1};
Line(21) = {13, 5};
Line(22) = {14, 6};
Line(23) = {10, 2};
Line(24) = {12, 4};
Line(25) = {11, 3};
Line(26) = {15, 17};
Line(27) = {17, 7};
Line(28) = {16, 18};
Line(29) = {18, 8};
Line(30) = {18, 17};
Line Loop(31) = {1, 2, 3, 4, 5, 6};
Plane Surface(32) = {31};
Line Loop(33) = {5, -24, -15, 25};
Plane Surface(34) = {33};
Line Loop(35) = {4, -25, -14, 23};
Plane Surface(36) = {35};
Line Loop(37) = {23, -3, -22, 13};
Plane Surface(38) = {37};
Line Loop(39) = {22, -2, -21, 12};
Plane Surface(40) = {39};
Line Loop(41) = {1, -21, -10, 20};
Plane Surface(42) = {41};
Line Loop(43) = {6, -20, -16, 24};
Plane Surface(44) = {43};
Line Loop(45) = {13, 14, 15, 16, 10, 12};
Plane Surface(46) = {45};
Line Loop(47) = {19, -12, 17, 18};
Plane Surface(48) = {47};
Line Loop(49) = {22, -9, -29, -28, 19};
Plane Surface(50) = {49};
Line Loop(51) = {9, -2, 7, 8};
Plane Surface(52) = {51};
Line Loop(53) = {7, -27, -26, -17, 21};
Plane Surface(54) = {53};
Line Loop(55) = {8, -29, 30, 27};
Plane Surface(56) = {55};
Line Loop(57) = {28, 30, -26, 18};
Plane Surface(58) = {57};
Surface Loop(59) = {46, 38, 36, 32, 42, 44, 34, 40};
Volume(60) = {59};
Surface Loop(61) = {54, 52, 50, 56, 58, 48, 40};
Volume(62) = {61};

Physical Surface("MOVING")    = {32};
Physical Surface("NOSLIP")    = {38, 42, 54, 52, 50};
Physical Surface("OUTFLOWL")  = {44};
Physical Surface("OUTFLOWR")  = {36};
Physical Surface("OUTFLOWB")  = {34};
Physical Surface("INLET1")    = {56};
Physical Surface("INLET2")    = {58};
Physical Surface("PNULL")     = {48, 46};
Physical Volume("FLUID")      = {62, 60};

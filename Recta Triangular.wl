(*Suposición:los eslabones están hecho de MDF*)\[Rho] = 
  500;(*densidad MDF=500[kg/m^3]*)
g = -9.81;(*Aceleración gravitatoria*)

(*Datos del eslabón 1 o Cilindro*)
d1 = 15;(*Valor de d1/Altura del cilindro*)
RadCil = 2.5;(*Radio Cilindro*)
m[1] = \[Rho] (\[Pi] (RadCil/100)^2*(d1/100));(*masa del cilindro*)

(*Datos del eslabón 2*)
e2 = 5;(*Espesor eslabón 2*)
g2 = 5;(*Ancho eslabón 2*)
l2 = 25;(*Largo eslabón 2*)
m[2] = \[Rho] (e2/100) (g2/100) (l2/100);(*masa eslabón 2*)

(*Datos del eslabón 3*)
e3 = 5;(*Espesor eslabón 3*)
g3 = 5;(*Ancho eslabón 3*)
l3 = 20;(*Largo eslabón 3*)
m[3] = \[Rho] (e3/100) (g3/100) (l3/100);(*masa eslabón 3*)


(*------------------------------------------------------------------------------------------------------------------*)
(*------------------------------------------------------------------------------------------------------------------*)
(*------------------------------------------------------------------------------------------------------------------*)

(*Construcción del plano de trabajo*)
NA = {-15, 30, 25};
NB = {-15, 20, 15};
NC = {15, 20, 15};
ND = {15, 30, 25};
plano = Line[{NA, NB, NC, ND, NA}];

(*Construcción de la trayectoria en línea recta*)
Linea = Line[{NA, NC}];
P0 = NA;
Pf = NC;

(*Longitud de la recta*)
qf = Sqrt[(Pf[[1]] - P0[[1]])^2 + (Pf[[2]] - P0[[2]])^2] // N;
(*Tiempo de operación*)
tf = 5;

(*Perfil de trayectoria triangular para q*)
q = Piecewise[{{(2*qf/tf^2)*t^2, t < tf/2}, 
               {(-2*qf/tf^2)*t^2 + (4*qf/tf)*t - qf, t >= tf/2}}];

(*Parametrizando la recta*)
x = P0[[1]] + q/qf (Pf[[1]] - P0[[1]]);
y = P0[[2]] + q/qf (Pf[[2]] - P0[[2]]);
z = P0[[3]] + q/qf (Pf[[3]] - P0[[3]]);

(*Generar puntos de la línea recta*)
puntosLinea = Table[{x, y, z}, {t, 0, tf, tf/50}];

(*Línea de la trayectoria con etiquetas*)
lineaTrayectoria = Line[puntosLinea];

(*Añadir etiquetas*)
etiquetas = {Text["NA", NA, {2, 0}], Text["NB", NB, {2, 0}], 
   Text["NC", NC, {-2, 0}], Text["ND", ND, {-2, 0}]};

Graphics3D[{plano, Red, lineaTrayectoria, etiquetas}]


(*------------------------------------------------------------------------------------------------------------------*)
(*------------------------------------------------------------------------------------------------------------------*)
(*------------------------------------------------------------------------------------------------------------------*)

(*Longitud eslabón 1*)
l1 = d1;

(*l2 y l3 conocidos desde el inicio*)
(*l4,l5,l6,\[Alpha],\[Beta] t \[Gamma] vienen del análisis de la \
cinemática inversa*)

l4 = z - l1;
l5 = Sqrt[x^2 + y^2];
l6 = Sqrt[l4^2 + l5^2];

\[Theta][1] = ArcTan[x, y];

\[Alpha] = ArcTan[l5, l4];
\[Beta] = 
  ArcTan[l2^2 + l6^2 - l3^2, 
   Sqrt[4 l2^2 l6^2 - (l2^2 + l6^2 - l3^2)^2]];

\[Theta][2] = \[Alpha] + \[Beta];

\[Gamma] = 
  ArcTan[l2^2 + l3^2 - l6^2, 
   Sqrt[4 l2^2 l3^2 - (l2^2 + l3^2 - l6^2)^2]];
\[Theta][3] = \[Gamma] + \[Pi];

(*Derivadas articulares*)
\[Theta]p[1] = D[\[Theta][1], t];
\[Theta]p[2] = D[\[Theta][2], t];
\[Theta]p[3] = D[\[Theta][3], t];
\[Theta]pp[1] = D[\[Theta]p[1], t];
\[Theta]pp[2] = D[\[Theta]p[2], t];
\[Theta]pp[3] = D[\[Theta]p[3], t];

(******************)

(*Modelado Directo*)
T01 = ({{Cos[\[Theta][1]], 0, Sin[\[Theta][1]], 0}, {Sin[\[Theta][1]],
      0, -Cos[\[Theta][1]], 0}, {0, 1, 0, d1}, {0, 0, 0, 1}});
T12 = ({{Cos[\[Theta][2]], -Sin[\[Theta][2]], 0, 
     l2 Cos[\[Theta][2]]}, {Sin[\[Theta][2]], Cos[\[Theta][2]], 0, 
     l2 Sin[\[Theta][2]]}, {0, 0, 1, 0}, {0, 0, 0, 1}});
T23 = ({{Cos[\[Theta][3]], -Sin[\[Theta][3]], 0, 
     l3 Cos[\[Theta][3]]}, {Sin[\[Theta][3]], Cos[\[Theta][3]], 0, 
     l3 Sin[\[Theta][3]]}, {0, 0, 1, 0}, {0, 0, 0, 1}});
H = ({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}});

(*Parámetros para tunear*)
ejes = Axes -> True;
aspecto = AspectRatio -> 1/1;
rango = PlotRange -> {{-50, 50}, {-50, 50}, {-20, 50}};
letreros = AxesLabel -> {"x", "y", "z"};
origen = AxesOrigin -> {0, 0, 0};
pov = ViewPoint -> {0.7, -1.5, 1};
(***************************)
(*Construcción eslabón 1*)
(*Separación del centro de la base respecto al marco de trabajo*)
Sep = 0;
(*Centro de la base*)
CenBase = {0, Sep, 0};
(*Centro de la tapa*)
CenTapa = {0, Sep, d1};
eslabon1 = Cylinder[{CenBase, CenTapa}, RadCil];
(*Construcción eslabón 2*)
(*Puntos en coordenadas homogéneas*)
V1 = {0, 1/2 g2, 1/2 e2, 1};
V2 = {0, -(1/2) g2, 1/2 e2, 1};
V3 = {0, -(1/2) g2, -(1/2) e2, 1};
V4 = {0, 1/2 g2, -(1/2) e2, 1};
V5 = {-l2, 1/2 g2, 1/2 e2, 1};
V6 = {-l2, -(1/2) g2, 1/2 e2, 1};
V7 = {-l2, -(1/2) g2, -(1/2) e2, 1};
V8 = {-l2, 1/2 g2, -(1/2) e2, 1};
B1 = H . T01 . T12 . V1;
B2 = H . T01 . T12 . V2;
B3 = H . T01 . T12 . V3;
B4 = H . T01 . T12 . V4;
B5 = H . T01 . T12 . V5;
B6 = H . T01 . T12 . V6;
B7 = H . T01 . T12 . V7;
B8 = H . T01 . T12 . V8;
cara1 = Polygon[{B1, B2, B3, B4, B1}];
cara2 = Polygon[{B5, B1, B2, B6, B5}];
cara3 = Polygon[{B6, B2, B3, B7, B6}];
cara4 = Polygon[{B8, B4, B3, B7, B8}];
cara5 = Polygon[{B5, B1, B4, B8, B5}];
cara6 = Polygon[{B5, B6, B7, B8, B5}];
eslabon2 = {cara1, cara2, cara3, cara4, cara5, cara6};

(*Construcción eslabón 3*)
N1 = {0, 1/2 g3, 1/2 e3, 1};
N2 = {0, -(1/2) g3, 1/2 e3, 1};
N3 = {0, -(1/2) g3, -(1/2) e3, 1};
N4 = {0, 1/2 g3, -(1/2) e3, 1};
N5 = {-l3, 1/2 g3, 1/2 e3, 1};
N6 = {-l3, -(1/2) g3, 1/2 e3, 1};
N7 = {-l3, -(1/2) g3, -(1/2) e3, 1};
N8 = {-l3, 1/2 g3, -(1/2) e3, 1};
M1 = H . T01 . T12 . T23 . N1;
M2 = H . T01 . T12 . T23 . N2;
M3 = H . T01 . T12 . T23 . N3;
M4 = H . T01 . T12 . T23 . N4;
M5 = H . T01 . T12 . T23 . N5;
M6 = H . T01 . T12 . T23 . N6;
M7 = H . T01 . T12 . T23 . N7;
M8 = H . T01 . T12 . T23 . N8;
kara1 = Polygon[{M1, M2, M3, M4, M1}];
kara2 = Polygon[{M5, M1, M2, M6, M5}];
kara3 = Polygon[{M6, M2, M3, M7, M6}];
kara4 = Polygon[{M8, M4, M3, M7, M8}];
kara5 = Polygon[{M5, M1, M4, M8, M5}];
kara6 = Polygon[{M5, M6, M7, M8, M5}];
eslabon3 = {kara1, kara2, kara3, kara4, kara5, kara6};

imagen = 
  Graphics3D[{Pink, eslabon1, Yellow, eslabon2, Cyan, eslabon3, Red, 
    plano, Purple, Linea}, ejes, rango, letreros, origen, pov];
Tabla = Table[imagen, {t, 0, tf, .01}];

ListAnimate[Tabla]

(*Gráficas de los pares reactivos*)
\[Tau]m1 = 
  1.5/100 (1/
      12 (6 g Cos[\[Theta][2]] (l3 Cos[\[Theta][3]] m[3] + 
          l2 (m[2] + 2 m[3])) - 
       6 l3 m[3] Sin[\[Theta][3]] (g Sin[\[Theta][2]] + 
          l2 (\[Theta]p[2] (2 \[Theta]p[1] + \[Theta]p[2]) - 
             2 (\[Theta]p[1] + \[Theta]p[2]) \[Theta]p[3] - \[Theta]p[
               3]^2)) + 
       4 l2^2 m[2] \[Theta]pp[
         1] + (3 (l1^2 + 2 RadCil^2) m[1] + g2^2 m[2]) \[Theta]pp[1] +
        g3^2 m[3] \[Theta]pp[1] + 12 l2^2 m[3] \[Theta]pp[1] + 
       4 l3^2 m[3] \[Theta]pp[1] + g2^2 m[2] \[Theta]pp[2] - 
       2 l2^2 m[2] \[Theta]pp[2] + g3^2 m[3] \[Theta]pp[2] + 
       4 l3^2 m[3] \[Theta]pp[2] + 
       6 l2 l3 Cos[\[Theta][3]] m[
         3] (2 \[Theta]pp[1] + \[Theta]pp[2] - \[Theta]pp[3]) + 
       g3^2 m[3] \[Theta]pp[3] - 2 l3^2 m[3] \[Theta]pp[3]));
\[Tau]m2 = 
  1.5/100 (1/
      12 (-6 g l2 Cos[\[Theta][2]] m[2] + 
       6 g l3 Cos[\[Theta][2] + \[Theta][3]] m[3] + 
       6 l2 l3 m[3] Sin[\[Theta][3]] \[Theta]p[1]^2 + 
       g2^2 m[2] \[Theta]pp[1] - 2 l2^2 m[2] \[Theta]pp[1] + 
       g3^2 m[3] \[Theta]pp[1] + 4 l3^2 m[3] \[Theta]pp[1] + 
       6 l2 l3 Cos[\[Theta][3]] m[3] \[Theta]pp[1] + 
       g2^2 m[2] \[Theta]pp[2] + 4 l2^2 m[2] \[Theta]pp[2] + 
       g3^2 m[3] \[Theta]pp[2] + 4 l3^2 m[3] \[Theta]pp[2] + 
       g3^2 m[3] \[Theta]pp[3] - 2 l3^2 m[3] \[Theta]pp[3]));
\[Tau]m3 = 
  1.5/100 (1/12 m[
      3] (-2 l3 (3 g Cos[\[Theta][2] + \[Theta][3]] + 
          3 l2 Sin[\[Theta][3]] \[Theta]p[1]^2 + 
          3 l2 Cos[\[Theta][3]] \[Theta]pp[1] + 
          l3 (\[Theta]pp[1] + \[Theta]pp[2] - 2 \[Theta]pp[3])) + 
       g3^2 (\[Theta]pp[1] + \[Theta]pp[2] + \[Theta]pp[3])));



(*Gráfica de la evolución del perfil Triangular*)
Plot[q, {t, 0, tf}, AxesLabel -> {"t", "Perfil Triangular"}, 
 PlotStyle -> Blue, AxesStyle -> Directive[Thick, Black]]

(*Cálculo de la velocidad del perfil Triangular*)
velocidadPerfil = D[q, t];

(*Gráfica de la velocidad del perfil Triangular*)
Plot[velocidadPerfil, {t, 0, tf}, 
 AxesLabel -> {"t", "Velocidad del Perfil Triangular"}, 
 PlotStyle -> Red, AxesStyle -> Directive[Thick, Black]]

 (*Cálculo de la aceleración del perfil Triangular*)
aceleracionPerfil = D[velocidadPerfil, t];

 (*Gráfica de Aceleración del Perfil Triangular*)
Plot[aceleracionPerfil, {t, 0, tf}, 
 AxesLabel -> {"t", "Aceleración del Perfil Triangular"}, 
 PlotStyle -> Green, AxesStyle -> Directive[Thick, Black]]


(*Gráfica del desplazamiento articular de cada articulación*)
Plot[{\[Theta][1] (180/\[Pi]), \[Theta][2] (180/\[Pi]), \[Theta][
    3] (180/\[Pi])}, {t, 0, tf}, 
 AxesLabel -> {"t", "Ángulo (grados)"}, 
 PlotLegends -> {"Articulación 1", "Articulación 2", 
   "Articulación 3"}, PlotStyle -> {Blue, Red, Green}, 
 AxesStyle -> Directive[Thick, Black]]

(*Gráfica de las derivadas articulares (velocidad) de cada \
articulación*)
Plot[{\[Theta]p[1] (180/\[Pi]), \[Theta]p[2] (180/\[Pi]), \[Theta]p[
    3] (180/\[Pi])}, {t, 0, tf}, 
 AxesLabel -> {"t", "Velocidad (grados/segundo)"}, 
 PlotLegends -> {"Articulación 1", "Articulación 2", 
   "Articulación 3"}, PlotStyle -> {Blue, Red, Green}, 
 AxesStyle -> Directive[Thick, Black]]

(*Gráfica de las derivadas articulares (aceleración) de cada \
articulación*)
Plot[{\[Theta]pp[1] (180/\[Pi]), \[Theta]pp[
    2] (180/\[Pi]), \[Theta]pp[3] (180/\[Pi])}, {t, 0, tf}, 
 AxesLabel -> {"t", "Aceleración (grados/segundo^2)"}, 
 PlotLegends -> {"Articulación 1", "Articulación 2", 
   "Articulación 3"}, PlotStyle -> {Blue, Red, Green}, 
 AxesStyle -> Directive[Thick, Black]]

(*Gráfica de los pares reactivos*)
Plot[{\[Tau]m1, \[Tau]m2, \[Tau]m3}, {t, 0, tf}, 
 AxesLabel -> {"t", "Par Reactivo"}, 
 PlotLegends -> {"Articulación 1", "Articulación 2", 
   "Articulación 3"}, PlotStyle -> {Blue, Red, Green}, 
 AxesStyle -> Directive[Thick, Black]]
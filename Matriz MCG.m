%% Dinámica del Robot usando Euler-Lagrange
%% Hallando matriz de inercia (M)
clc,
syms m1 m2 m3 m4 m5 m6 lc1 lc2 lc3 lc4 lc5 lc6 l1 l2 l3 l4 l5 l6 g0_1 g0_2 g0_3...
q1 q2 q3 q4 q5 q6 dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6...
Ixx1 Ixy1 Ixz1 Iyx1 Iyy1 Iyz1 Izx1 Izy1 Izz1 Ixx2 Ixy2 Ixz2 Iyx2 Iyy2 Iyz2 Izx2 Izy2 Izz2...
Ixx3 Ixy3 Ixz3 Iyx3 Iyy3 Iyz3 Izx3 Izy3 Izz3 Ixx4 Ixy4 Ixz4 Iyx4 Iyy4 Iyz4 Izx4 Izy4 Izz4 ...
Ixx5 Ixy5 Ixz5 Iyx5 Iyy5 Iyz5 Izx5 Izy5 Izz5 Ixx6 Ixy6 Ixz6 Iyx6 Iyy6 Iyz6 Izx6 Izy6 Izz6 g

% 1.- Hallar las posiciones de los centros de masa

clc,
% xc1 = lc1*cos(q1); 
% yc1 = lc1*sin(q1); 
% zc1 = 0;
% xc2 = l1*cos(q1)+lc2*cos(q1+q2);
% yc2 = l1*sin(q1)+lc2*sin(q1+q2); 
% zc2 = 0;
xc1=(0.205/2)*cos(q1);
xc1 = vpa(xc1,4)

yc1=(0.205/2)*sin(q1);
yc1 = vpa(yc1,4)

xc2=0.205*cos(q1)+0.05*cos(q1+q2);
xc2 = vpa(xc2,4)

yc2=0.205*sin(q1)+0.05*sin(q1+q2);
yc2 = vpa(yc2,4)

xc3=0.205*cos(q1)+0.1*cos(q1+q2)-0.038*cos(q1+q2+q3);
xc3 = vpa(xc3,4)

yc3=0.205*sin(q1)+0.1*sin(q1+q2)-0.038*sin(q1+q2+q3);
yc3 = vpa(yc3,4)

xc4=0.205*cos(q1)+0.1*cos(q1+q2)-(0.076+q4)*cos(q1+q2+q3);
xc4 = vpa(xc4,4)

yc4=0.205*sin(q1)+0.1*sin(q1+q2)-(0.076+q4)*sin(q1+q2+q3);
yc4 = vpa(yc4,4)

zc1 = 0;
zc2 = 0;
zc3 = 0;
zc4 = 0;
% 2.- Jacobianos de velocidad lineal
% OJO: VER SI SON DE REVOLUCIÓN O PRISMÁTICA
n = 4; % # Articulaciones

Jv1 = vpa(Jv_CdM(xc1,yc1,zc1,1,n),4) % R
Jv2 = vpa(Jv_CdM(xc2,yc2,zc2,2,n),4); % R
Jv2 = vpa(Jv2(1:3,1:4),4)
Jv3 = vpa(Jv_CdM(xc3,yc3,zc3,3,n),4) % R
Jv4 = [0 0 0 0;
       0 0 0 0;
       1 1 1 1]
% 3.- Jacobianos angulares 
% ver "Dyn Manipuladores I" diapo25 para mas detalles
z0 = [0;0;1];
z1 = [0;0;1];
z2 = [0;0;1];
z3 = [0;0;1];

Jw1 = [z0 zeros(3,3)]
Jw2 = [z0 z1 zeros(3,2)] 
Jw3 = [z0 z1 z2 zeros(3,1)]
Jw4 = zeros(3,4)

% 4.- Productos de Jacobianos de velocidad lineal
 
Jprod_v_1 = simplify(Jv1.'*Jv1);
Jprod_v_1 = vpa(Jprod_v_1,4)

Jprod_v_2 = simplify(Jv2.'*Jv2);
Jprod_v_2 = vpa(Jprod_v_2,4)

Jprod_v_3 = simplify(Jv3.'*Jv3);
Jprod_v_3 = vpa(Jprod_v_3,4)

Jprod_v_4 = Jv4.'*Jv4
% 5.- Productos de Jacobianos de velocidad angular con R

% OJO: LOS R DEBEN ESTAR CON RESPECTO AL ORIGEN
T01 = DHtoMatrix(0,q1,0.205,0) 
T12 = DHtoMatrix(0,q2-(pi/2),-0.1,0)
T23 = DHtoMatrix(0,q3,0,pi/2);
T23(1,2)=-0.0007963*sin(q3);
T23(2,2)=0.0007963*cos(q3);
T23(3,3)=0.0007963;
T23 = vpa(simplify(T23),2)
T34 = DHtoMatrix(q4+0.117,0,0,0)

T02 = vpa(simplify(T01*T12),4)

T03 = vpa(simplify(T02*T23),4)
T04 = vpa(simplify(T03*T34),4)

R01 = vpa(simplify(T01(1:3,1:3)),4)
R02 = vpa(simplify(T02(1:3,1:3)),4)
R03 = vpa(simplify(T03(1:3,1:3)),4)
R04 = vpa(simplify(T04(1:3,1:3)),4)

Jw1R1 = Jw1.'*R01 
Jw2R2 = Jw2.'*R02 
Jw3R3 = vpa(Jw3.'*R03 ,4)
Jw4R4 = Jw4.'*R04 

% 6.- Productos del Tensor de Inercia de cada CdM
I1 = [Ixx1 Ixy1 Ixz1; Iyx1 Iyy1 Iyz1; Izx1 Izy1 Izz1];
I2 = [Ixx2 Ixy2 Ixz2; Iyx2 Iyy2 Iyz2; Izx2 Izy2 Izz2];
I3 = [Ixx3 Ixy3 Ixz3; Iyx3 Iyy3 Iyz3; Izx3 Izy3 Izz3];
I4 = [Ixx4 Ixy4 Ixz4; Iyx4 Iyy4 Iyz4; Izx4 Izy4 Izz4];

I1=[6.31E+06 1.21E-08 -2.84E-08;
    1.21E-08 6.01E+06 3.18E-07;
    -2.84E-08 3.18E-07 5.35E+05]
I2=[3.24E+05 2.49E-08 0;
    2.49E-08 1.41E+06 -1.56E-08;
    0 -1.56E-08 1.59E+06]
I3=[5.15E+05 -2.49E-09 1.921;
    -2.49E-09 3.90E+05 1.07E-08;
    1.921 1.07E-08 2.11E+05]
I4=[79043.028 -3.55E-10 7.76E-06;
    -3.55E-10 79142.642 4.44E-09;
    7.76E-06 4.44E-09 8424.281]

pt1 = Jw1R1*I1*R01.'*Jw1
pt2 = Jw2R2*I2*R02.'*Jw2
pt3 = vpa(Jw3R3*I3*R03.'*Jw3,4)
pt4 = Jw4R4*I4*R04.'*Jw4

% 7.- Reemplazando términos (ver diapo51 para más detalles)

% básicamente es M = (mi*Jvi'*Jvi) + (Jwi'*Ri*Ii*Ri'*Jwi)
% es decir, multiplicar todo lo del paso 4 pero multiplicado con sus masas respectivas
ParteLineal = m1*Jprod_v_1 + m2*Jprod_v_2 +  m3*Jprod_v_3 + m4*Jprod_v_4;
ParteLineal = vpa(ParteLineal,4)
% luego sumar todo lo del paso 6, sin multiplicar nada a parte
ParteAngular = pt1 + pt2 + pt3 + pt4;
ParteAngular = vpa(ParteAngular,4)

% Por ultimo, sumar lineal y angular

M = ParteLineal + ParteAngular;
M = vpa(M,4)

%% Hallando Coriolis (C)
clc,
qvector = [q1 q2 q3 q4 q5 q6];
dqvector = [dq1 dq2 dq3 dq4 dq5 dq6];
n = 4; %articulaciones

% Como son n gdl, se tendrán n^3 c (c111,c112,...cnnn)
% la cantidad de pares repetidos son (n^2 - n)*n/2
% con ello, tendríamos que hallar n^3 - (n^2 - n)*n/2 términos.
% Pero esa operación se resumió (al parecer, omitió) en una sola linea, hallado de frente c_ij para la matriz C 


c = 0; 
for i=1:n
    for j=1:n
        for k=1:n
            pt1 = diff( M(i,j) , qvector(1,k) );
            pt2 = diff( M(i,k) , qvector(1,j) );
            pt3 = diff( M(j,k) , qvector(1,i) );
            
            c = c + ( (1/2)*( pt1 + pt2 - pt3 ) ) * dqvector(1,k);
        end
        C(i,j) = c;
        c = 0;
    end
end
C = vpa(simplify(C),4)
%% Hallando Gravedad (G)
clc,
%g0 = [g0_1 g0_2 g0_3]; % Casi siempre será [0 -9.81 0]
g0 = [0; -g; 0];
m = [m1 m2 m3 m4 m5 m6];
Jv_vector =  [Jv1 Jv2 Jv3 Jv4]; % Editar en caso cuantos Jv tengas
% dimensiones de G: n x 1 
% depende cuantas articulaciones tengas para poder hallar dicha matriz

g1 = -Jv1(:,1).'*m1*g0 - Jv2(:,1).'*m2*g0 - Jv3(:,1).'*m3*g0 - Jv4(:,1).'*m4*g0 ;
g2 = -Jv1(:,2).'*m1*g0 - Jv2(:,2).'*m2*g0 - Jv3(:,2).'*m3*g0 - Jv4(:,2).'*m4*g0 ;
g3 = -Jv1(:,3).'*m1*g0 - Jv2(:,3).'*m2*g0 - Jv3(:,3).'*m3*g0 - Jv4(:,3).'*m4*g0 ;
g4 = -Jv1(:,4).'*m1*g0 - Jv2(:,4).'*m2*g0 - Jv3(:,4).'*m3*g0 - Jv4(:,4).'*m4*g0 ;

% Darse cuenta de como va aumentando (:,1)
g1 = vpa(simplify(g1),4)
g2 = vpa(simplify(g2),4)
g3 = vpa(simplify(g3),4)
g4 = vpa(simplify(g4),4)



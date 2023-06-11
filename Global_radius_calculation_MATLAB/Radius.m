clc
clear all
close all

% -------------- coordinates of 3 points -----------------------
% i th point: (x(i),y(i),z(i))
x=[8.4540 22.6728 37.8185];
y=[36.7525 31.3499 27.5957];
z=[424.9901 426.9355 425.3460];
% -------------- Coplanar constraint -----------------------
gm=[991    992    993    1;
    x(1)   y(1)   z(1)   1;
    x(2)   y(2)   z(2)   1;
    x(3)   y(3)   z(3)   1];
abc0=det(gm)*inv(gm); % Adjoint matrix of gm
A1=abc0(1,1);
B1=abc0(2,1);
C1=abc0(3,1);
D1=abc0(4,1);
% -------------- Distance Equality Constraint -------------------
A2=2*(x(2)-x(1));
B2=2*(y(2)-y(1));
C2=2*(z(2)-z(1));
D2=x(1)^2+y(1)^2+z(1)^2-x(2)^2-y(2)^2-z(2)^2;
A3=2*(x(3)-x(1));
B3=2*(y(3)-y(1));
C3=2*(z(3)-z(1));
D3=x(1)^2+y(1)^2+z(1)^2-x(3)^2-y(3)^2-z(3)^2;
% ----------- Calculate coordinates of the circle center-------------
ABC=[A1 B1 C1;
     A2 B2 C2;
     A3 B3 C3];
XYZ=-inv(ABC)*[D1;D2;D3];
% -------------- Calculate radius -----------------------
m=(x(1)-XYZ(1))^2+(y(1)-XYZ(2))^2+(z(1)-XYZ(3))^2;
R=sqrt(m);
disp(R); % export radius



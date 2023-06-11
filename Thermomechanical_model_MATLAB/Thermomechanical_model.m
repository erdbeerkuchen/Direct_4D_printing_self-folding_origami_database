clc
clear all
close all

% -------------- input parameter ---------------------------------------
% --- Geometry ---
t1_ratio = 10/20;             % thickness ratio of the first layer
theta = [-45 45]*pi/180;      % angle of each layer
alpha = 90;                   % critical angle of the diamond
diag = 16;                    % long diagonal of the panel element
lgedg = sqrt(diag^2/(tand(alpha/2)));
l = lgedg/sqrt(2);             % initial length [mm]
w = lgedg/sqrt(2);             % initial width  [mm]

% --- Material ---
% E1,E2,G12: modulus [Pa];   v12,v21: Poisson's Ratio;    e12_T:  thermal strain 

E1 = 8.26; v12 = 0.5; v21 = 0.4; E2 = E1*v21/v12; G12 = 0.79; e12_T = [-0.05787; 0.00681; 0]; % T=65
% E1 = 5.82; v12 = 0.5; v21 = 0.4; E2 = E1*v21/v12; G12 = 0.78; e12_T = [-0.10680; 0.02568; 0]; % T=75
% E1 = 4.87; v12 = 0.5; v21 = 0.4; E2 = E1*v21/v12; G12 = 0.67; e12_T = [-0.14539; 0.02997; 0]; % T=85  

% ----------------------------------------------------------------------

t_total = 2.0;            % total thickness [mm]
t1 = t_total*t1_ratio;
t2 = t_total-t1;
t = [t1 t2];               % thickness of each layer
n = length(t);              % number of layers

%% Classical laminate theory to calculate the strain and curvature

Q = [ E1/(1-v12*v21)       E1*v21/(1-v12*v21)       0;
     E1*v21/(1-v12*v21)     E2/(1-v12*v21)          0;
         0                        0                G12];
     
%%calculate the z-coordanates 
m = sum(t)/2;
M = zeros(1,n+1);
for loopx=1:n
    M(loopx+1) = M(loopx)+t(loopx);
end
z = M-m;     

%%thermal resultant force and moment
N_T = 0;
M_T = 0;

A = 0;
B = 0;
D = 0;

% e12_T = 0.5*e12_T; %check model

%%assemble 
for loopx=1:n
    % transformation matrix
    T(:,:,loopx) = [cos(theta(loopx))^2  sin(theta(loopx))^2  2*sin(theta(loopx))*cos(theta(loopx));
         sin(theta(loopx))^2  cos(theta(loopx))^2  -2*sin(theta(loopx))*cos(theta(loopx));
         -sin(theta(loopx))*cos(theta(loopx))  sin(theta(loopx))*cos(theta(loopx)) cos(theta(loopx))^2-sin(theta(loopx))^2;];
    Qbar(:,:,loopx) = inv(T(:,:,loopx))*Q*(inv(T(:,:,loopx)))';
    
    exy_T = (T(:,:,loopx)')*e12_T;
    
    N_T = N_T+Qbar(:,:,loopx)*exy_T*(z(loopx+1)-z(loopx));
    M_T = M_T+1/2*Qbar(:,:,loopx)*exy_T*(z(loopx+1)^2-z(loopx)^2);
    
    A = A+Qbar(:,:,loopx)*(z(loopx+1)-z(loopx));
    B = B+1/2*Qbar(:,:,loopx)*(z(loopx+1)^2-z(loopx)^2);
    D = D+1/3*Qbar(:,:,loopx)*(z(loopx+1)^3-z(loopx)^3);
    
end

%calculation

ABD = [A B;B D];
invABD = inv([A B;B D]);

%initialize full-screen figure
figure('Units', 'normalized', 'Position', [0 0 1 1]);

for loopN = 1:1:1
    %%increment the loading
    N_Tloc = N_T*loopN;
    M_Tloc = M_T*loopN;

    %%strains at this increment
    sol = invABD*[N_Tloc;M_Tloc];
    e0 = sol(1:3,1);            % strain at mid-plane
    k  = sol(4:6,1);            % curvature

    %%displacements from strains
    px = e0(1);                 %exx
    qy = px;                    %eyy
    rxy = k(3)/2;                 %kxy

    %symmetric
    A = [l l 0; l/2 l/2 1; l -l 0];%l and w are equal at this point
    b = [-rxy*(l)^2; -rxy*(l/2)^2; 0];

    solu = inv(A)*b;
    rx = solu(1);
    ry = solu(2);
    r0 = solu(3);

    %self written function for aspect ratio
    coeff = [px qy rxy rx ry r0];
    [xu, yv, ww, xx, yy] = trimesh_func(alpha, diag, coeff);
    calb_ww = max(max(ww));

    %prepare half and quarter surfaces
    ww_half1 = ww;
    ww_half2 = ww;
    ww_half3 = ww;
    ww_half4 = ww;
    ww_q1 = ww;
    ww_q2 = ww;
    ww_q3 = ww;
    ww_q4 = ww;
    for loopx = 1:length(xx)
        for loopy = 1:length(yy)
            if loopx <= loopy
                ww_half1(loopx, loopy) = nan;
            elseif loopx >= loopy
                ww_half2(loopx, loopy) = nan;
            end
            if (loopx+loopy) >= length(xx)
                ww_half3(loopx, loopy) = nan;
            elseif (loopx+loopy) <= length(xx)
                ww_half4(loopx, loopy) = nan;
            end
            if loopx <= loopy || (loopx+loopy) >= length(xx)
                ww_q1(loopx, loopy) = nan;
            end
            if loopx <= loopy || (loopx+loopy) <= length(xx)
                ww_q2(loopx, loopy) = nan;
            end
            if loopx >= loopy || (loopx+loopy) <= length(xx)
                ww_q3(loopx, loopy) = nan;
            end
            if loopx >= loopy || (loopx+loopy) >= length(xx)
                ww_q4(loopx, loopy) = nan;
            end
        end
    end

    %%plot surface
    %yoshimura geometric relationship
    [xid, yid] = find(ww==max(max(ww)));
    xid = xid(1);yid = yid(1);
    nd1 = [xu(2,2) yv(2,2) ww(2,2)];
    nd2 = [xu(xid,yid) yv(xid,yid) ww(xid,yid)];
    nd3 = [xu(500,500) yv(500,500) ww(500,500)];
    nd4 = [xu(yid,xid) yv(yid,xid) ww(yid,xid)];

    lgedg = norm(nd2 - nd4);%b
    stedg = norm(nd1 - nd3);%d
    sdlgth = norm(nd1 - nd2);%a


    pause(0.001);%pause and delete plots to demonstrate animation

    %yoshimura pattern
    theta_rot = atand(sqrt(4*sdlgth^2-lgedg^2-stedg^2)/lgedg)*4;
    %radius of formed cylindrical shell
    r_tot = lgedg/(2*tand(theta_rot/2)) + calb_ww;

    if loopN == 0 %infinite radius when flat
        %radius of formed cylindrical shell
        r_tot = lgedg/2/sqrt(3)*400 + calb_ww;%400 is empirical value
        theta_rot = atand(lgedg/2/(r_tot - calb_ww))*2;
    end

    %surface plotting
%     if loopN > 0
%         delete(s0);delete(s1);delete(s2);
%         delete(s3);delete(s4);delete(s5);
%         delete(s6);delete(s7);delete(s8);
%         delete(s9);delete(s10);delete(s11);
%         delete(s12);delete(s13);delete(s14);
%         delete(s15);delete(s16);delete(s17);
%         delete(s18);delete(s19);delete(s20);
%         delete(s21);delete(s22);delete(s23);
%         delete(s24);
%     end

    colormap("copper");
    clim("manual");
    clim([-2.24 2.24]);

    %layer#1 3 surfaces
    s0 = surf(xu, yv, ww, ww);hold on;
    s1 = surf(xu, yv, ww, ww);
    s2 = surf(xu, yv, ww, ww);

    direction = [1 1 0];
    origin = [xu(251,251) yv(251,251) r_tot];

    rotate(s0, direction, theta_rot, origin);
    rotate(s1, direction, theta_rot*2, origin);

    %layer#2 2 surfaces and 2 half surfaces
    xincre = abs(xu(1,501) - xu(501,1));
    yincre = abs(yv(1,501) - yv(501,1));
    s3 = surf(xu+xincre/2, yv+yincre/2, ww, ww);hold on;
    s4 = surf(xu+xincre/2, yv+yincre/2, ww, ww);
    s5 = surf(xu+xincre/2, yv+yincre/2, ww_half2, ww);
    s9 = surf(xu+xincre/2, yv+yincre/2, ww_half1, ww);

    direction = [1 1 0];
    origin = [xu(251,251)+xincre/2 yv(251,251)+yincre/2 r_tot];

    rotate(s3, direction, theta_rot*0.5, origin);
    rotate(s4, direction, theta_rot*1.5, origin);
    rotate(s5, direction, theta_rot*2.5, origin);
    rotate(s9, direction, -theta_rot*0.5, origin);

    %layer#3 3 surfaces
    s6 = surf(xu+xincre, yv+yincre, ww, ww);hold on;
    s7 = surf(xu+xincre, yv+yincre, ww, ww);
    s8 = surf(xu+xincre, yv+yincre, ww, ww);

    direction = [1 1 0];
    origin = [xu(251,251)+xincre yv(251,251)+yincre r_tot];

    rotate(s6, direction, theta_rot, origin);
    rotate(s7, direction, theta_rot*2, origin);

    %layer#4 2 half surfaces and 2 quarter surfaces
    s10 = surf(xu+xincre*1.5, yv+yincre*1.5, ww, ww);
    s11 = surf(xu+xincre*1.5, yv+yincre*1.5, ww, ww);
    s12 = surf(xu+xincre*1.5, yv+yincre*1.5, ww_half2, ww);
    s13 = surf(xu+xincre*1.5, yv+yincre*1.5, ww_half1, ww);

    direction = [1 1 0];
    origin = [xu(251,251)+xincre*1.5 yv(251,251)+yincre*1.5 r_tot];

    rotate(s10, direction, theta_rot*0.5, origin);
    rotate(s11, direction, theta_rot*1.5, origin);
    rotate(s12, direction, theta_rot*2.5, origin);
    rotate(s13, direction, -theta_rot*0.5, origin);

    %layer#5 3 surfaces
    s18 = surf(xu+xincre*2, yv+yincre*2, ww_half3, ww);hold on;
    s19 = surf(xu+xincre*2, yv+yincre*2, ww_half3, ww);
    s20 = surf(xu+xincre*2, yv+yincre*2, ww_half3, ww);

    direction = [1 1 0];
    origin = [xu(251,251)+xincre*2.0 yv(251,251)+yincre*2.0 r_tot];

    rotate(s18, direction, theta_rot, origin);
    rotate(s19, direction, theta_rot*2, origin);

    %layer#6 2 half surfaces and 2 quarter surfaces
    s21 = surf(xu-xincre*0.5, yv-yincre*0.5, ww_half1, ww);
    s22 = surf(xu-xincre*0.5, yv-yincre*0.5, ww, ww);
    s23 = surf(xu-xincre*0.5, yv-yincre*0.5, ww, ww);
    s24 = surf(xu-xincre*0.5, yv-yincre*0.5, ww_half2, ww);

    direction = [1 1 0];
    origin = [xu(251,251)-xincre*0.5 yv(251,251)-yincre*0.5 r_tot];

    rotate(s21, direction, -theta_rot*0.5, origin);
    rotate(s22, direction, theta_rot*0.5, origin);
    rotate(s23, direction, theta_rot*1.5, origin);
    rotate(s24, direction, theta_rot*2.5, origin);

    
    %layer#7 2 half surfaces and 2 quarter surfaces
    s14 = surf(xu-xincre, yv-yincre, ww_half4, ww);
    s15 = surf(xu-xincre, yv-yincre, ww_half4, ww);
    s16 = surf(xu-xincre, yv-yincre, ww_half4, ww);

    direction = [1 1 0];
    origin = [xu(251,251)-xincre yv(251,251)-yincre r_tot];

    rotate(s15, direction, theta_rot*1.0, origin);
    rotate(s16, direction, theta_rot*2.0, origin);

    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    shading interp;
    axis([-80 80 -40 120 -40 80]);

    pause(0.001);%pause and delete plots to demonstrate animation
end




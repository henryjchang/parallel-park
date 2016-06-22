yalmip('clear')
clear all

% Model Parameters (Do not touch)
lr = 0.13;
dr = 0.06;
lf = 0.13;
df = 0.07;
w = 0.095;
Lf = lf+df;
Lr = lr+dr;
L_a = lf; %m (distance to front axis)
L_b = lr; %m (distance to back axis)
length = Lf+Lr; % m
width = w*2; % m

% Horizon
N1 = 5;
N2 = 7;
% steps
T1 = 75;
T2 = 40;
dt = 0.2;

% Initial Pose
x0 = [-1; 0.6; 0; 0];

% Target pose
xT = [0;w+.06;0;0];
x_ref = xT(1);
y_ref = xT(2);
psi_ref = xT(3);
v_ref = xT(4);

% Intermediate pose
xI = [xT(1) + 0.4; xT(2) + 0.225 ; 0; 0];
%xI = [-0.5; 0.6 ; 0; 0];

% Parking spot definition, state and input constraints
xl = xT(1) - (Lr + .07); 
xr = xT(1) + Lf + .1;
yt = xT(2) + w;

% Parking spot size
spot_l = xr - xl;
spot_w = yt;

zmax = [20;20;4*pi;10];
zmin = [-20;0;-4*pi;-10];
umax = [0.6;1.5*dt];
x_min = zmin(1);
y_min = zmin(2);
psi_min = zmin(3);
v_min = zmin(4);
x_max = zmax(1);
y_max = zmax(2);
psi_max = zmax(3);
v_max = zmax(4);

nx = 4; % number of states
nu = 2; % number of inputs

%% First maneuver

%Model Constraints
a_max = 1.5;
a_min = -1.3;
v_max = 2;
v_min = -1*v_max;
d_f_max = pi/6;
d_f_min = -1*d_f_max;


%MPC data
N = 7; %horizon
dt = 0.2; %time step
T = 35; %number of steps

u = sdpvar(repmat(nu,1,N),repmat(1,1,N)); % (nu x 1) x N
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1)); % (nx x 1) x (N+1)
%bta = zeros(N,1); % Nx1

objective = 0;


constraints = [];
for k = 1:N 
    %Rajamani p. 24
    constraints = [constraints,
                    %bta(k)       == atan(L_b/(L_a+L_b)*tan(u{k}(1))),
                    x{k+1}(1)    == x{k}(1) + dt*(x{k}(4)*cos(x{k}(3) + atan(L_b/(L_a+L_b)*tan(u{k}(1))))), % x
                    x{k+1}(2)    == x{k}(2) + dt*(x{k}(4)*sin(x{k}(3) + atan(L_b/(L_a+L_b)*tan(u{k}(1))))), % y
                    x{k+1}(3)    == x{k}(3) + dt*(x{k}(4)*cos(atan(L_b/(L_a+L_b)*tan(u{k}(1))))/(L_a+L_b)*tan(u{k}(1))), % psi
                    x{k+1}(4)    == x{k}(4) + dt*u{k}(2), % v
                    %v_min <= x{k}(4) <= v_max,
                    d_f_min <= u{k}(1) <= d_f_max, 
                    a_min <= u{k}(2) <= a_max
                  ];
    %objective = objective + (x{k}(1) - x_ref)^2 + (x{k}(2) - y_ref)^2 + ...
    %                   (x{k}(3) - psi_ref)^2 + (x{k}(4) - v_ref)^2;          
end

objective = objective + (x{N+1}(1) - xI(1))^2 + (x{N+1}(2) - xI(2))^2 + ...
                        (x{N+1}(3) - xI(3))^2 + (x{N+1}(4) - xI(4))^2;


parameters_in = {x{1}};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','ipopt'),parameters_in,solutions_out);


oldx = x0;


%use these for animation purposes
x_vals = [];
y_vals = [];
psi_vals = [];


axis([-1.5 1 -.2 1.5]);
hold on
for t1 = 1:T
    [solutions,diagnostics] = controller{{oldx}};
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    iteration = t1
    
    X = solutions{2};
    oldx = X(:,2);
    
    x_vals = [x_vals, oldx(1)];
    y_vals = [y_vals, oldx(2)];
    psi_vals = [psi_vals, oldx(3)];
    
    x_vertices = [x_vals(t1) + Lf*cos(psi_vals(t1)) - width/2*sin(psi_vals(t1));
                  x_vals(t1) + Lf*cos(psi_vals(t1)) + width/2*sin(psi_vals(t1));
                  x_vals(t1) - Lr*cos(psi_vals(t1)) + width/2*sin(psi_vals(t1));
                  x_vals(t1) - Lr*cos(psi_vals(t1)) - width/2*sin(psi_vals(t1));
                  x_vals(t1) + Lf*cos(psi_vals(t1)) - width/2*sin(psi_vals(t1));
                  ];
    y_vertices = [y_vals(t1) + Lf*sin(psi_vals(t1)) + width/2*cos(psi_vals(t1));
                  y_vals(t1) + Lf*sin(psi_vals(t1)) - width/2*cos(psi_vals(t1));
                  y_vals(t1) - Lr*sin(psi_vals(t1)) - width/2*cos(psi_vals(t1));
                  y_vals(t1) - Lr*sin(psi_vals(t1)) + width/2*cos(psi_vals(t1));
                  y_vals(t1) + Lf*sin(psi_vals(t1)) + width/2*cos(psi_vals(t1));
                  ];
           
    plot(x_vertices, y_vertices, 'k');
    pause(dt);
    hold on;
    
    if norm(oldx - xI, 2) <= 0.2
        break;
    end
end
%% Second maneuver

u2 = sdpvar(repmat(nu,1,N2),repmat(1,1,N2)); % (nu x 1) x N2
x2 = sdpvar(repmat(nx,1,N2+1),repmat(1,1,N2+1)); % (nx x 1) x (N2+1)
d = binvar(repmat(16,1,N2),repmat(1,1,N2));


constraints2 = [];
for k = 1:N2
    %Rajamani p. 24
    constraints2 = [constraints2,
                    x2{k+1}(1)    == x2{k}(1) + dt*(x2{k}(4)*cos(x2{k}(3) + u2{k}(1))), % x
                    x2{k+1}(2)    == x2{k}(2) + dt*(x2{k}(4)*sin(x2{k}(3) + u2{k}(1))), % y
                    x2{k+1}(3)    == x2{k}(3) + dt*(x2{k}(4)/lr*sin(u2{k}(1))), % psi
                    x2{k+1}(4)    == x2{k}(4) + dt*u2{k}(2), % v
                    -umax <= u2{k} <= umax,
                    zmin <= x2{k} <= zmax,
                  ];

%     x_vert =     [x{k}(1) + length/2*cos(x{k}(3)) - width/2*sin(x{k}(3));
%                   x{k}(1) + length/2*cos(x{k}(3)) + width/2*sin(x{k}(3));
%                   x{k}(1) - length/2*cos(x{k}(3)) + width/2*sin(x{k}(3));
%                   x{k}(1) - length/2*cos(x{k}(3)) - width/2*sin(x{k}(3));
%                   ];
%     y_vert =     [x{k}(2) + length/2*sin(x{k}(3)) + width/2*cos(x{k}(3));
%                   x{k}(2) + length/2*sin(x{k}(3)) - width/2*cos(x{k}(3));
%                   x{k}(2) - length/2*sin(x{k}(3)) - width/2*cos(x{k}(3));
%                   x{k}(2) - length/2*sin(x{k}(3)) + width/2*cos(x{k}(3));
%                   ];

% Parking spot obstacle constraints              
% x >= xl OR y >= yt
% x <= xr OR y >= yt
% evauluated at all four corners of car

    constraints2 = [constraints2, d{k}(1)*(x2{k}(1) + Lf*cos(x2{k}(3)) - width/2*sin(x2{k}(3)) - xl) + ...
        d{k}(2)*(x2{k}(2) + Lf*sin(x2{k}(3)) + width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(1) + d{k}(2) == 1];
    constraints2 = [constraints2, d{k}(3)*(x2{k}(1) + Lf*cos(x2{k}(3)) + width/2*sin(x2{k}(3)) - xl) + ...
        d{k}(4)*(x2{k}(2) + Lf*sin(x2{k}(3)) - width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(3) + d{k}(4) == 1];
    constraints2 = [constraints2, d{k}(5)*(x2{k}(1) - Lr*cos(x2{k}(3)) + width/2*sin(x2{k}(3)) - xl) + ... 
        d{k}(6)*(x2{k}(2) - Lr*sin(x2{k}(3)) - width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(5) + d{k}(6) == 1];
    constraints2 = [constraints2, d{k}(7)*(x2{k}(1) - Lr*cos(x2{k}(3)) - width/2*sin(x2{k}(3)) - xl) + ...
        d{k}(8)*(x2{k}(2) - Lr*sin(x2{k}(3)) + width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(7) + d{k}(8) == 1];
    
    
    constraints2 = [constraints2, d{k}(9)*(-(x2{k}(1) + Lf*cos(x2{k}(3)) - width/2*sin(x2{k}(3))) + xr) + ...
        d{k}(10)*(x2{k}(2) + Lf*sin(x2{k}(3)) + width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(9) + d{k}(10) == 1];
    constraints2 = [constraints2, d{k}(11)*(-(x2{k}(1) + Lf*cos(x2{k}(3)) + width/2*sin(x2{k}(3))) + xr) + ... 
        d{k}(12)*(x2{k}(2) + Lf*sin(x2{k}(3)) - width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(11) + d{k}(12) == 1];
    constraints2 = [constraints2, d{k}(13)*(-(x2{k}(1) - Lr*cos(x2{k}(3)) + width/2*sin(x2{k}(3))) + xr) + ...
        d{k}(14)*(x2{k}(2) - Lr*sin(x2{k}(3)) - width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(13) + d{k}(14) == 1];
    constraints2 = [constraints2, d{k}(15)*(-(x2{k}(1) - Lr*cos(x2{k}(3)) - width/2*sin(x2{k}(3))) + xr) + ...
        d{k}(16)*(x2{k}(2) - Lr*sin(x2{k}(3)) + width/2*cos(x2{k}(3)) - yt) >= 0,
        d{k}(15) + d{k}(16) == 1];
    
    constraints2 = [constraints2, y_min <= x2{k}(2) + Lf*sin(x2{k}(3)) + width/2*cos(x2{k}(3)) <= y_max,
        y_min <= x2{k}(2) + Lf*sin(x2{k}(3)) - width/2*cos(x2{k}(3)) <= y_max,
        y_min <= x2{k}(2) - Lr*sin(x2{k}(3)) - width/2*cos(x2{k}(3)) <= y_max,
        y_min <= x2{k}(2) - Lr*sin(x2{k}(3)) + width/2*cos(x2{k}(3)) <= y_max];
    constraints2 = [constraints2, x_min <= x2{k}(1) + Lf*cos(x2{k}(3)) - width/2*sin(x2{k}(3)) <= x_max,
        x_min <= x2{k}(1) + Lf*cos(x2{k}(3)) + width/2*sin(x2{k}(3)) <= x_max,
        x_min <= x2{k}(1) - Lr*cos(x2{k}(3)) + width/2*sin(x2{k}(3)) <= x_max,
        x_min <= x2{k}(1) - Lr*cos(x2{k}(3)) - width/2*sin(x2{k}(3)) <= x_max];
   
end

%constraints = [constraints, x2{N+1} == [x_ref;y_ref;psi_ref;v_ref]];
%Think about normalizing?

objectiveT = 100*((x2{N2+1}(1) - xT(1))^2 + (x2{N2+1}(2) - xT(2))^2 + ...
                        (x2{N2+1}(3) - xT(3))^2) + (x2{N2+1}(4) - xT(4))^2;

parameters_in2 = {x2{1}};
solutions_out2 = {[u2{:}], [x2{:}]};


controllerT = optimizer(constraints2, objectiveT,sdpsettings('solver','ipopt'),parameters_in2,solutions_out2);



for t2 = 1:T2
    [solutions,diagnostics] = controllerT{{oldx}};
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    iteration = t2
    
    X = solutions{2};
    oldx = X(:,2);
    
    x_vals = [x_vals, oldx(1)];
    y_vals = [y_vals, oldx(2)];
    psi_vals = [psi_vals, oldx(3)];
    
    x_vertices = [x_vals(t2) + Lf*cos(psi_vals(t2)) - width/2*sin(psi_vals(t2));
                  x_vals(t2) + Lf*cos(psi_vals(t2)) + width/2*sin(psi_vals(t2));
                  x_vals(t2) - Lr*cos(psi_vals(t2)) + width/2*sin(psi_vals(t2));
                  x_vals(t2) - Lr*cos(psi_vals(t2)) - width/2*sin(psi_vals(t2));
                  x_vals(t2) + Lf*cos(psi_vals(t2)) - width/2*sin(psi_vals(t2));
                  ];
    y_vertices = [y_vals(t2) + Lf*sin(psi_vals(t2)) + width/2*cos(psi_vals(t2));
                  y_vals(t2) + Lf*sin(psi_vals(t2)) - width/2*cos(psi_vals(t2));
                  y_vals(t2) - Lr*sin(psi_vals(t2)) - width/2*cos(psi_vals(t2));
                  y_vals(t2) - Lr*sin(psi_vals(t2)) + width/2*cos(psi_vals(t2));
                  y_vals(t2) + Lf*sin(psi_vals(t2)) + width/2*cos(psi_vals(t2));
                  ];
           
    plot(x_vertices, y_vertices, 'k');
    pause(dt);
    hold on;
    
    if norm(oldx - xT,2) <= 0.1
        break;
    end
end


%% Animation
x_vertices_ref = [x_ref + length/2*cos(psi_ref) - width/2*sin(psi_ref);
                  x_ref + length/2*cos(psi_ref) + width/2*sin(psi_ref);
                  x_ref - length/2*cos(psi_ref) + width/2*sin(psi_ref);
                  x_ref - length/2*cos(psi_ref) - width/2*sin(psi_ref);
                  x_ref + length/2*cos(psi_ref) - width/2*sin(psi_ref);
                  ];
y_vertices_ref = [y_ref + length/2*sin(psi_ref) + width/2*cos(psi_ref);
                  y_ref + length/2*sin(psi_ref) - width/2*cos(psi_ref);
                  y_ref - length/2*sin(psi_ref) - width/2*cos(psi_ref);
                  y_ref - length/2*sin(psi_ref) + width/2*cos(psi_ref);
                  y_ref + length/2*sin(psi_ref) + width/2*cos(psi_ref);
                  ]; 
figure()              
plot(x_vertices_ref, y_vertices_ref, 'r');
axis([-1 1 -.2 1]);
hold on;
x_boundary = [-1000, xl, xl, xr,...
                xr, 1000];
y_boundary = [yt, yt, 0, ...
                0, yt, yt];
plot(x_boundary, y_boundary);

for j = 1:T1+T2
    x_vertices = [x_vals(j) + Lf*cos(psi_vals(j)) - width/2*sin(psi_vals(j));
                  x_vals(j) + Lf*cos(psi_vals(j)) + width/2*sin(psi_vals(j));
                  x_vals(j) - Lr*cos(psi_vals(j)) + width/2*sin(psi_vals(j));
                  x_vals(j) - Lr*cos(psi_vals(j)) - width/2*sin(psi_vals(j));
                  x_vals(j) + Lf*cos(psi_vals(j)) - width/2*sin(psi_vals(j));
                  ];
    y_vertices = [y_vals(j) + Lf*sin(psi_vals(j)) + width/2*cos(psi_vals(j));
                  y_vals(j) + Lf*sin(psi_vals(j)) - width/2*cos(psi_vals(j));
                  y_vals(j) - Lr*sin(psi_vals(j)) - width/2*cos(psi_vals(j));
                  y_vals(j) - Lr*sin(psi_vals(j)) + width/2*cos(psi_vals(j));
                  y_vals(j) + Lf*sin(psi_vals(j)) + width/2*cos(psi_vals(j));
                  ];
           
    plot(x_vertices, y_vertices, 'k');
    pause(dt);
    hold on;
end

yalmip('clear')
clear all


lr = 0.13;
dr = 0.06;
lf = 0.13;
df = 0.07;
w = 0.095;
L = lr+lf;
Lf = lf+df;
Lr = lr+dr;




%Initial
z0 = [0.7;0.3;0;0];
zT = [0;w+.03;0;0];
obs = [-(Lr+0.01);Lf+.07;2*w+0.01];
zmax = [20;20;4*pi;10];
zmin = [-20;0;-4*pi;-10];
T = 21;
dt = 0.2;
umax = [0.6;1.5*dt];

x_min = zmin(1);
y_min = zmin(2);
psi_min = zmin(3);
v_min = zmin(4);

x_max = zmax(1);
y_max = zmax(2);
psi_max = zmax(3);
v_max = zmax(4);


%Targets
x_ref = zT(1);
y_ref = zT(2);
psi_ref = zT(3);
v_ref = zT(4);


nx = 4; % number of states
nu = 2; % number of inputs

%Model Parameters
L_a = lf; %m (distance to front axis)
L_b = lr; %m (distance to back axis)
% 
%Rectangular Approximation
length = Lf+Lr; % m
width = w*2; % m

% % Parking spot size
% % spot_l = 0.6; % m
% % spot_w = 0.22; % m

xl = obs(1);
xr = obs(2);
yt = obs(3);
spot_l = xr - xl;
spot_w = yt;


u = sdpvar(repmat(nu,1,T),repmat(1,1,T)); % (nu x 1) x T
x = sdpvar(repmat(nx,1,T+1),repmat(1,1,T+1)); % (nx x 1) x (T+1)
%bta = zeros(T,1); % Tx1
d = binvar(repmat(16,1,T),repmat(1,1,T));

objective = 0;
constraints = [];
for k = 1:T 
    %Rajamani p. 24
    constraints = [constraints,
                    %bta(k)       == atan(L_b/(L_a+L_b)*tan(u{k}(1))),
                    x{k+1}(1)    == x{k}(1) + dt*(x{k}(4)*cos(x{k}(3) + u{k}(1))), % x
                    x{k+1}(2)    == x{k}(2) + dt*(x{k}(4)*sin(x{k}(3) + u{k}(1))), % y
                    x{k+1}(3)    == x{k}(3) + dt*(x{k}(4)/lr*sin(u{k}(1))), % psi
                    x{k+1}(4)    == x{k}(4) + dt*u{k}(2), % v
                    %v_min <= x{k}(4) <= v_max,
                    -umax <= u{k} <= umax,
                    zmin <= x{k} <= zmax,
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

    constraints = [constraints, d{k}(1)*(x{k}(1) + Lf*cos(x{k}(3)) - width/2*sin(x{k}(3)) - xl) + ...
        d{k}(2)*(x{k}(2) + Lf*sin(x{k}(3)) + width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(1) + d{k}(2) == 1];
    constraints = [constraints, d{k}(3)*(x{k}(1) + Lf*cos(x{k}(3)) + width/2*sin(x{k}(3)) - xl) + ...
        d{k}(4)*(x{k}(2) + Lf*sin(x{k}(3)) - width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(3) + d{k}(4) == 1];
    constraints = [constraints, d{k}(5)*(x{k}(1) - Lr*cos(x{k}(3)) + width/2*sin(x{k}(3)) - xl) + ... 
        d{k}(6)*(x{k}(2) - Lr*sin(x{k}(3)) - width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(5) + d{k}(6) == 1];
    constraints = [constraints, d{k}(7)*(x{k}(1) - Lr*cos(x{k}(3)) - width/2*sin(x{k}(3)) - xl) + ...
        d{k}(8)*(x{k}(2) - Lr*sin(x{k}(3)) + width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(7) + d{k}(8) == 1];
    
    
    constraints = [constraints, d{k}(9)*(-(x{k}(1) + Lf*cos(x{k}(3)) - width/2*sin(x{k}(3))) + xr) + ...
        d{k}(10)*(x{k}(2) + Lf*sin(x{k}(3)) + width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(9) + d{k}(10) == 1];
    constraints = [constraints, d{k}(11)*(-(x{k}(1) + Lf*cos(x{k}(3)) + width/2*sin(x{k}(3))) + xr) + ... 
        d{k}(12)*(x{k}(2) + Lf*sin(x{k}(3)) - width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(11) + d{k}(12) == 1];
    constraints = [constraints, d{k}(13)*(-(x{k}(1) - Lr*cos(x{k}(3)) + width/2*sin(x{k}(3))) + xr) + ...
        d{k}(14)*(x{k}(2) - Lr*sin(x{k}(3)) - width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(13) + d{k}(14) == 1];
    constraints = [constraints, d{k}(15)*(-(x{k}(1) - Lr*cos(x{k}(3)) - width/2*sin(x{k}(3))) + xr) + ...
        d{k}(16)*(x{k}(2) - Lr*sin(x{k}(3)) + width/2*cos(x{k}(3)) - yt) >= 0,
        d{k}(15) + d{k}(16) == 1];
    
    constraints = [constraints, y_min <= x{k}(2) + Lf*sin(x{k}(3)) + width/2*cos(x{k}(3)) <= y_max,
        y_min <= x{k}(2) + Lf*sin(x{k}(3)) - width/2*cos(x{k}(3)) <= y_max,
        y_min <= x{k}(2) - Lr*sin(x{k}(3)) - width/2*cos(x{k}(3)) <= y_max,
        y_min <= x{k}(2) - Lr*sin(x{k}(3)) + width/2*cos(x{k}(3)) <= y_max];
    constraints = [constraints, x_min <= x{k}(1) + Lf*cos(x{k}(3)) - width/2*sin(x{k}(3)) <= x_max,
        x_min <= x{k}(1) + Lf*cos(x{k}(3)) + width/2*sin(x{k}(3)) <= x_max,
        x_min <= x{k}(1) - Lr*cos(x{k}(3)) + width/2*sin(x{k}(3)) <= x_max,
        x_min <= x{k}(1) - Lr*cos(x{k}(3)) - width/2*sin(x{k}(3)) <= x_max];
    
%     constraints = [constraints, x_ref - 5 <= x{k}(1) <= x_ref + 5, ...
%                         y_ref - spot_w/2 <= x{k}(2) <= y_ref + 3];
    %objective = objective + (x{T}(1) - x_ref)^2 + (x{T}(2) - y_ref)^2;
    
end

constraints = [constraints, x{T+1} == [x_ref;y_ref;psi_ref;v_ref]];
%Think about normalizing?
objective = objective + 100*((x{T+1}(1) - x_ref)^2 + (x{T+1}(2) - y_ref)^2 + ...
                        (x{T+1}(3) - psi_ref)^2) + ...
                        (x{T+1}(4) - v_ref)^2;
% objective = 0;
% for t = T-4:T
%     objective = objective + (x{t}(3) - zT(3))^2;
% end


parameters_in = {x{1}};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective,sdpsettings('solver','ipopt'),parameters_in,solutions_out);


oldx = z0;

%use these for animation purposes
x_vals = [];
y_vals = [];
psi_vals = [];

[solutions,diagnostics] = controller{{oldx}};    
if diagnostics == 1
        error('The problem is infeasible');
end
    
X = solutions{2};
x_vals = X(1,:);
y_vals = X(2,:);
psi_vals = X(3,:);




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

for j = 1:T
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
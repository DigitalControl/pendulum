%
% Control the inverted pendulum
%
%close all;
clear all;
pkg load control;

% Whether to create plots
plotAll = false;

%
% Pendulum model, the longer rod
%

% Sample rate / control frequency (Hz)
f = 300;
T = 1/f;
Maxpos = 0.25;              % Max carriage travel +- 0.25 m
Maxangle = 0.175;           % Max rod angle -- 10 deg
Maxvoltage = 20;            % Max motor voltage, V
pstart = 0.005;             % Carriage position starting limit, m
astart = 1*pi/180;          % Angle starting limit, rad
g = 9.81;                   % m/s^2     Gravitational constant

% SYSTEM PARAMETERS
% Measured Mechanical Parameters
d1 = 0.323;                % m         Length of pendulum 1 (long)
d2 = 0.079;                 % m         Length of pendulum 2 (short)
%mp1 = 0.0208;              % kg        Mass of pendulum 1
mp1 = 0.0318;
%mp2 = 0.0050;              % kg        Mass of pendulum 2
mp2 = 0.0085;
%m=0.8; % Makes state estimates better for some reason
m = 0.3163;                 % kg        Mass of carriage
rd = 0.0254/2;              % m         Drive pulley radius
md = 0.0375;                % kg        Mass of drive pulley (cylinder)
%mc1 = 0.0036;              % kg        Mass of clamp 1*
%mc2 = 0.0036;              % kg        Mass of clamp 2*
mc1 = 0.0085;
mc2 = mc1;

% *Clamp Dimensions
%  Rectangular 0.0254 x 0.01143 m
%  The pivot shaft is 0.00714 m from the end

% Motor Parameters (Data Sheet)
Im = 43e-7;                 % kg m^2/rad    Rotor moment of inertia
R = 4.09;                   % ohms          Resistance
kt = 0.0351;                % Nm/A          Torque constant
ke = 0.0351;                % Vs/rad        Back emf constant

% Derived Mechanical Parameters

                            % kg m^2/rad    Moment of inertia, clamp 1
%Ic1 = mc1*(0.01143^2 + 0.0254^2)/12 + mc1*(0.0127-0.00714)^2;
Ic1 = mc1*(0.0098^2 + 0.0379^2)/12;
Ic2 = Ic1;                  % kg m^2/rad    Moment of inertia, clamp 2
Id = md*(rd^2)/2;           % kg m^2/rad    Moment of inertia, drive pulley
Imd = Im + Id;              % kg m^2/rad    Moment of inertia, combined

J1 = Ic1 + mp1*(d1^2)/3;    % Total moment of inertia, pendulum 1 (long)
J2 = Ic2 + mp2*(d2^2)/3;    % Total moment of inertia, pendulum 2 (short)
Jd = Im + Id;               % Total moment of inertia, motor drive
Mc = m + mc1 + mc2;         % Total carriage mass

% Friction Test Data
%   Carriage Slope = 19 deg;  Terminal Velocity xdotss = 0.312 m/s; From
%       twincarriage.m; formula b = m g sin(theta)/xdotss
%   Pendulum 1 (long) Exponent a1 = 0.0756 1/s;  From longfit.m
%   Pendulum 2 (short) Exponent a2 = 0.2922 1/s; From shortfit.m
%        formula b = 2 a J

%alpha = 19;
alpha = 12.2;
%xdotss = 0.312;
xdotss = 0.4852;
%a1 = 0.0756;
%a2 = 0.2922;
a1 = 0.0185;
a2 = 0.012;
                            % Ns/m    Viscous friction of carriage system
b = (Mc + mp1 + mp2)*g*sin(alpha*pi/180)/xdotss;
b1 = 2*a1*J1;               % Nms/rad    Viscous friction of pendulum 1 (rotational)
b2 = 2*a2*J2;               % Nms/rad    Viscous friction of pendulum 2 (rotational)

scale = [rd*2*pi/4096  2*pi/4096 -0.05/250];

%
% System modeling
%
% Combining motor and two pendulum equations
%
l1 = d1/2;
l2 = d2/2;
Km = kt;
bc = 1; % not given, so guess

% Denominators p1, p2
p1 = (R*rd^2*(Mc*l1^2*l2^2*mp1*mp2+J1*Mc*l2^2*mp2+J1*l2^2*mp1*mp2+J2*Mc*l1^2*mp1+J2*l1^2*mp1*mp2+J1*J2*Mc+J1*J2*mp1+J1*J2*mp2));
p2 = (Mc*l1^2*l2^2*mp1*mp2+J1*Mc*l2^2*mp2+J1*l2^2*mp1*mp2+J2*Mc*l1^2*mp1+J2*l1^2*mp1*mp2+J1*J2*Mc+J1*J2*mp1+J1*J2*mp2);

A = [
    % xd
    0 1 0 0 0 0;
    % xdd
    0 -(R*bc*l1^2*l2^2*mp1*mp2*rd^2+J1*R*bc*l2^2*mp2*rd^2+J2*R*bc*l1^2*mp1*rd^2+Km^2*l1^2*l2^2*mp1*mp2+J1*J2*R*bc*rd^2+J1*Km^2*l2^2*mp2+J2*Km^2*l1^2*mp1+J1*J2*Km^2)/p1 -(-R*g*l1^2*l2^2*mp1^2*mp2*rd^2-J2*R*g*l1^2*mp1^2*rd^2)/p1 0 -(-R*g*l1^2*l2^2*mp1*mp2^2*rd^2-J1*R*g*l2^2*mp2^2*rd^2)/p1 0;
    % theta1d
    0 0 0 1 0 0;
    % theta1dd
    0 -(R*bc*l2^2*mp2*rd^2+J2*R*bc*rd^2+Km^2*l2^2*mp2+J2*Km^2)*l1*mp1/p1 -(-Mc*R*g*l2^2*mp2*rd^2-R*g*l2^2*mp1*mp2*rd^2-J2*Mc*R*g*rd^2-J2*R*g*mp1*rd^2-J2*R*g*mp2*rd^2)*l1*mp1/p1 0 +g*l2^2*mp2^2*l1*mp1/p2 0;
    % theta2d
    0 0 0 0 0 1;
    % theta2dd
    0 +mp2*l2*(R*bc*l1^2*mp1*rd^2+J1*R*bc*rd^2+Km^2*l1^2*mp1+J1*Km^2)/p1 -mp2*l2*g*l1^2*mp1^2/p2 0 +mp2*l2*(-Mc*R*g*l1^2*mp1*rd^2-R*g*l1^2*mp1*mp2*rd^2-J1*Mc*R*g*rd^2-J1*R*g*mp1*rd^2-J1*R*g*mp2*rd^2)/p1 0;
    ];
B = [0;
    -(-Km*l1^2*l2^2*mp1*mp2*rd-J1*Km*l2^2*mp2*rd-J2*Km*l1^2*mp1*rd-J1*J2*Km*rd)/p1;
     0;
    -(-Km*l2^2*mp2*rd-J2*Km*rd)*l1*mp1/p1;
     0;
     mp2*l2*(-Km*l1^2*mp1*rd-J1*Km*rd)/p1;
     ];
C = [1 0 0 0 0 0;
     0 0 1 0 0 0];
D = [0;
     0];

states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta1_dot'};
inputs = {'v'};
outputs = {'x'; 'phi'};
sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

% Poles for our system
original_poles = eig(A);

%
% Discretize
%
sys_d = c2d(sys_ss, T, 'zoh');

% Check if controllable and observable
co = ctrb(sys_d);
ob = obsv(sys_d);

controllability = rank(co);
observability = rank(ob);

if controllability == 6
	disp('System is controllable. Yay!')
else
	disp('System is not controllable. This is bad.')
end
if observability == 6
	disp('System is observable. Yay!')
else
	disp('System is not observable. This is bad.')
end

% Original poles are
%
poles = eig(A);

%
% Create a controller with LQR and simulate it
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace
%
Q = C'*C;
Q(1,1) = 1000;
Q(3,3) = 100;
R = 0.1;

K = lqr(A,B,Q,R);

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:4;
r = 0.2*ones(size(t));
if plotAll
    figure;
    [y,t,x]=lsim(sys_cl,r,t);
    [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('Step Response with LQR Control');
end

%
% Correcting the cart position error. Now the cart actually ends up at 0.2
% meters as we were commanding it.
%
Cn = [1 0 0 0 0 0];
sys_ss = ss(A,B,Cn,0);
Nbar = rscale(sys_ss,K);

sys_cl = ss(Ac,Bc*Nbar,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

r =0.2*ones(size(t));
if plotAll
    figure;
    [y,t,x]=lsim(sys_cl,r,t);
    [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('Step Response with Precompensation and LQR Control');
end

%
% Designing observer/estimator
%
poles = eig(Ac);

%P = [-30 -31 -32 -33 -34 -35];
P = [-2 -3 -4 -5 -6 -7];
L = place(A',C',P)';

%
%
%
Ace = [(A-B*K) (B*K);
       zeros(size(A)) (A-L*C)];
Bce = [B*Nbar;
       zeros(size(B))];
Cce = [Cc zeros(size(Cc))];
Dce = [0;0];

states_est = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta1_dot' 'e1' 'e2' 'e3' 'e4' 'e5' 'e6'};
sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states_est,'inputname',inputs,'outputname',outputs);

r = 0.2*ones(size(t));
if plotAll
    figure;
    [y,t,x]=lsim(sys_est_cl,r,t);
    [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('Step Response with Observer-Based State-Feedback Control');
end

%
% Simulating the plant separately from our control and esimation, getting
% closer to what will actually run on the pendulum
%

% To get discretized
sys_est_only = ss(A, [B L], C, [D [0;0] [0;0]]);
sys_est_only_d = c2d(sys_est_only, T, 'zoh');

if plotAll
%if false
    % We have 4 state variables, and need the current and past values
    state = zeros(6, 2);
    eststate = zeros(6, 2);

    % The states we want to plot
    N = 4*f;
    output = zeros(N, 2);
    estoutput = zeros(2, 1);
    estoutputhistory = zeros(N, 2);

    % Constant input command
    r = 0.2;
    input = zeros(N);

    for i = 1:N
        % Our control law
        input(i) = r*Nbar - K*eststate(:,1);
        %input(i) = r*Nbar - K*(eststate(:,1) - [1;0;0;0]*r);
        %input(i) = r*Nbar - K*[state(1,1); eststate(2,1); state(3,1); eststate(4,1)];

        % Max out at +/- Maxvoltage
        input(i) = min(max(input(i),-Maxvoltage),Maxvoltage);

        % Estimation
        y = [state(1,1); state(3,1)]; % Use only the measured part of the state
        uvec = [input(i); y-estoutput];
        eststate(:,1) = sys_est_only_d.a*eststate(:,2) + sys_est_only_d.b*uvec;
        estoutput = sys_est_only_d.c*eststate(:,2) + sys_est_only_d.d*uvec;

        % Simulate actual system
        state(:,1) = sys_d.a*state(:,2) + sys_d.b*input(i);
        output(i,:) = sys_d.c*state(:,2) + sys_d.d*input(i);

        % Save the new states in the old state
        eststate(:,2) = eststate(:,1);
        state(:,2) = state(:,1);

        % Save so we can plot
        estoutputhistory(i,:) = estoutput';
    end

    t = 0:T:(size(estoutputhistory,1)-1)/f;
    %figure;
    %[AX,H1,H2] = plotyy(t,output(:,1),t,output(:,2),'plot');
    %set(get(AX(1),'Ylabel'),'String','cart position (m)');
    %set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    %title('actual state - without lsim');
    figure;
    [AX,H1,H2] = plotyy(t,estoutputhistory(:,1),t,estoutputhistory(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('State Estimates - without lsim');

    % Plot the input, make sure it's not > Maxvoltage
    figure;
    plot(t,estoutputhistory(:,1),'-r',
         t,estoutputhistory(:,2),'-b',
         t,input,'-g');
    title('Voltage Output');
end

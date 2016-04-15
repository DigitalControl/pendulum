pkg load control;

% Double Pendulum Parameters (Tentative:  There are two pendulum setups, each with different parameters.  I'm not sure which these go to.)
% This script is for balancing only the long rod.

% Run parameters
%f = input('Control Frequency (Hz) = ');
%Trun = input('Run Time (s) = ');
%f=130;
f=1000;
Trun=30;

kmax = round(f*Trun);
T = 1/f;
Maxpos = 0.25;              % Max carriage travel +- 0.25 m
Maxangle = 0.175;           % Max rod angle -- 10 deg
Maxvoltage = 20;            % Max motor voltage, V
pstart = 0.005;             % Carriage position starting limit, m
astart = 1*pi/180;          % Angle starting limit, rad

g = 9.81;                   % m/s^2     Gravitational constant

% SYSTEM PARAMETERS
% Measured Mechanical Parameters
d1 = -0.323;                 % m         Length of pendulum 1 (long)
d2 = 0.079;                 % m         Length of pendulum 2 (short)
%mp1 = 0.0208;              % kg        Mass of pendulum 1
mp1 = 0.0318;
%mp2 = 0.0050;              % kg        Mass of pendulum 2
mp2 = 0.0085;
%m=.5;
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
R = 4.09;                   % ohms            Resistance
kt = 0.0351;                % Nm/A            Torque constant
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


T = 1/f;

% The data above comes from the fweb wiki.

M=Mc;                       %mass of cart
m=mp1;                      %mass of pendulum 1
b=b;                        %friction
l=d1/2;                     %length of pendulum
I=J1;                       %inertia of pendulum
%q=(M+m)*(l+m*l^2)-(m*l)^2;
%num=[m*l,0];                                %numerator for transfer function
%den=[q,b*(l+m*l^2),-m*g*l*(M+m),-b*m*g*l];  %denominator for transfer function
%[A,B,C,D]=tf2ss(num,den)
%A,B,C,D matricies for the state space model
% x_vec is [x,  x_dot, theta, theta_dot]'
% See the web site:  https://www.library.cmu.edu/ctms/ctms/examples/pend/invpen.htm
A=[ 0   1                                                  0                                   0;
    0   ((-(I+m*l^2)*b)/(I*(M+m)+M*m*l^2))  ((m^2*g*l^2)/(I*(M+m)+M*m*l^2))     0;
    0   0                                   0                                   1;
    0   ((-m*l*b)/(I*(M+m)+M*m*l^2))        ((m*g*l*(M+m))/(I*(M+m)+M*m*l^2))   0];
B=[ 0;
    ((I+m*l^2)/(I*(M+m)+M*m*l^2));
    0;
    ((m*l)/(I*(M+m)+M*m*l^2))];
C=[ 1   0   0   0;
    0   0   1   0];
D=[ 0;
    0];
cont_sys = ss(A,B,C,D)
rank_ctrb = rank(ctrb(A,B))
original_poles = eig(A) %poles for our system
T=0.01
disc_sys = c2d(cont_sys, T)
H = tf(disc_sys)
%[num, den] = tf(disc_sys)

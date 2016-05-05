%
% Control the inverted pendulum
%
close all;
clear all;
pkg load control;

% Whether to create plots
plotAll = false;

%
% Pendulum model, the longer rod
%

% Sample rate / control frequency (Hz)
f = 100;
T = 1/f;
Maxpos = 0.25;              % Max carriage travel +- 0.25 m
Maxangle = 0.175;           % Max rod angle -- 10 deg
Maxvoltage = 20;            % Max motor voltage, V
pstart = 0.005;             % Carriage position starting limit, m
astart = 1*pi/180;          % Angle starting limit, rad
g = 9.81;                   % m/s^2     Gravitational constant

% SYSTEM PARAMETERS
% Measured Mechanical Parameters
d1 = -0.323;                % m         Length of pendulum 1 (long)
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
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
%
M = Mc;                     % Mass of cart
m = mp1;                    % Mass of pendulum 1
b = b;                      % Friction
I = J1;                     % Inertia of pendulum
g = 9.81;                   % m/s^2     Gravitational constant
l = d1/2;                   % Length of pendulum

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'u'};
outputs = {'x'; 'phi'};
sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

% Poles for our system
original_poles = eig(A);

%
% Discretize
%
sys_d = c2d(sys_ss, T, 'zoh');
%H = tf(sys_d);

% Check if controllable and observable
co = ctrb(sys_d);
ob = obsv(sys_d);

controllability = rank(co);
observability = rank(ob);

if controllability == 4
	disp('System is controllable. Yay!')
else
	disp('System is not controllable. This is bad.')
end
if observability == 4
	disp('System is observable. Yay!')
else
	disp('System is not observable. This is bad.')
end

% Original poles are:
%
%   0.00000 + 0.00000i
%  -4.44367 + 0.00000i
%  -0.04895 + 5.15645i
%  -0.04895 - 5.15645i
%
% which are way out of the unit circle.
%
poles = eig(A);

%
% Create a controller with LQR and simulate it
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace
%

% Try 1, cart flies away and the settling time is high (from example)
Q = C'*C;
R = 1;

% Try 2, cart doesn't fly away but the settling time is still high (from example)
Q = C'*C;
Q(1,1) = 5000;
Q(3,3) = 100;
R = 1;

% Try 3, better
Q = C'*C;
Q(1,1) = 2000;
Q(3,3) = 500;
R = 1;

K = lqr(A,B,Q,R);

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'r'};
outputs = {'x'; 'phi'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
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
Cn = [1 0 0 0];
sys_ss = ss(A,B,Cn,0);
Nbar = rscale(sys_ss,K);

sys_cl = ss(Ac,Bc*Nbar,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
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
% Our slowest pole's real value is at -1.7891. Let's try placing our
% estimator's poles at -20.
%
%  -10.0491 + 10.4058i
%  -10.0491 - 10.4058i
%   -1.7891 +  3.5523i
%   -1.7891 -  3.5523
%
poles = eig(Ac);

P = [-20 -21 -22 -23];
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

states = {'x' 'x_dot' 'phi' 'phi_dot' 'e1' 'e2' 'e3' 'e4'};
inputs = {'r'};
outputs = {'x'; 'phi'};

sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:5;
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
% Performing the simulation without using lsim
%
sys_est_d = c2d(sys_est_cl, T, 'zoh');

% This is slow and not really needed... just a proof of concept
if false
    % We have 8 state variables, and need the current and past values
    % Variables: x xdot phi phidot e1 e2 e3 e4
    state = zeros(8, 2);

    % The states we want to plot
    N = 2000;
    output = zeros(N, 2);

    % Constant input command of zero
    r = 0.2
    input = r*ones(N);
    %input = r*Nbar + K*state(1:4,2);

    for i = 1:N
        % We measure two of these states
        %state(1,2) = cartpos
        %state(3,2) = pendangle

        %input(i) = input(i) - K*state(1:4,2);
        %input(i) = L*(state(1:4,2) - state(5:8,2));

        state(:,1) = sys_est_d.a*state(:,2) + sys_est_d.b*input(i);

        % Save for plot
        output(i,:) = sys_est_d.c*state(:,2) + sys_est_d.d*input(i);

        % New values are now the old values
        state(:,2) = state(:,1);
    end

    figure;
    plot(output);
    t = 1:N;
    [AX,H1,H2] = plotyy(t,output(:,1),t,output(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('without lsim');
end

%
% Simulating the plant separately from our control and esimation, getting
% closer to what will actually run on the pendulum
%

% To get discretized: xbar(k+1) = (A-L*C)*xbar(k) + (B-L*D)*u + L*y
states = {'x' 'x_dot' 'phi' 'phi_dot'};
inputs = {'r'};
outputs = {'x'; 'phi'};
% Maybe B should be (B-L*D)*Nbar, but D is zero, so in this case it's the same
sys_est_only = ss(A-L*C,B*Nbar-L*D,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
sys_est_only_d = c2d(sys_est_only, T, 'zoh');
% Alt form
% See: http://cats-fs.rpi.edu/~wenj/ECSE444F01/discreteobserver.pdf
sys_est_alt = ss(A+L*(C-D*K)-B*K,L,-K,0);
sys_est_alt_d = c2d(sys_est_alt, T, 'zoh');

if false
    % We have 4 state variables, and need the current and past values
    state = zeros(4, 2);
    eststate = zeros(4, 2);

    % The states we want to plot
    N = 400;
    output = zeros(N, 2);
    estoutput = zeros(2, 1);
    estoutputhistory = zeros(N, 2);

    % Constant input command
    r = 0.2;
    input = zeros(N);

    for i = 1:N
        % Our control law
        input(i) = r*Nbar - K*state(:,1);
        %input(i) = r*Nbar - K*eststate(:,1);

        % Estimation
        % Same as actual, but with the additional L*y term at the end for
        % correction that we'll call F here
        y = [state(1,1) state(3,1)]'; % Use only the measured part of the state
        F = L*(y-estoutput);
        %F = L*C*state(:,1); % Also produces the same result using the entire state
        eststate(:,1) = sys_est_only_d.a*eststate(:,2) + sys_est_only_d.b*input(i) + F;
        estoutput = sys_est_only_d.c*eststate(:,2) + sys_est_only_d.d*input(i);

        %eststate(:,1) = sys_est_alt_d.a*eststate(:,2) + sys_est_alt_d.b*y + (B-L*D)*r*Nbar;
        %estoutput = [eststate(1,2) eststate(3,2)]';
        %estoutput = sys_est_alt_d.c*eststate(:,2) + r*Nbar;

        % Simulate actual system
        state(:,1) = sys_d.a*state(:,2) + sys_d.b*input(i);
        output(i,:) = sys_d.c*state(:,2) + sys_d.d*input(i);

        % Save the new states in the old state
        eststate(:,2) = eststate(:,1);
        state(:,2) = state(:,1);

        % Save so we can plot
        estoutputhistory(i,:) = estoutput';
    end

    figure;
    t = 1:N;
    [AX,H1,H2] = plotyy(t,output(:,1),t,output(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('actual state - without lsim');
    figure;
    [AX,H1,H2] = plotyy(t,estoutputhistory(:,1),t,estoutputhistory(:,2),'plot');
    set(get(AX(1),'Ylabel'),'String','cart position (m)');
    set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
    title('state estimates - without lsim');
end


%
% TF from cart position to PWM value
%
s = tf('s');
inductance = 0.01; % H, a complete guess, no idea... see if it works
H = (kt/rd)*(s*Jd+b1)/(s^2*Jd*inductance + s*(Jd*R+b1*inductance) + (b1*R+kt^2));
Hd = c2d(H, T, 'zoh');
% But, we really want the inverse, from F to V
% Note: this probably doesn't work... causes error if we c2d 1/H
Hd = 1/Hd;
[num, den] = tfdata(Hd);
num = cell2mat(num);
den = cell2mat(den);
b0 = num(1);
b1 = num(2);
b2 = num(3);
a0 = 0;
a1 = den(1);
a2 = den(2);

%x(1) =
%x_dot(1) =
%phi(1) =
%phi_dot(1) =

%e_x
%e_x_dot
%e_phi
%e_phi_dot


%H = tf(sys_d);

% Transfer function 'H' from input 'r' to output ...
%
%     6.312e-05 z^3 - 6.361e-05 z^2 - 6.212e-05 z + 6.262e-05
% x:  -------------------------------------------------------
%          z^4 - 3.976 z^3 + 5.929 z^2 - 3.93 z + 0.9766
%
%       -0.0001674 z^3 + 0.0001687 z^2 + 0.0001648 z - 0.0001661
% phi:  --------------------------------------------------------
%            z^4 - 3.976 z^3 + 5.929 z^2 - 3.93 z + 0.9766
%
%
% Converting these to the difference equations:
%
%    H_x(z) = X(z)/R(z)
% => X(z)(z^4 - 3.976 z^3 + 5.929 z^2 - 3.93 z + 0.9766) =
%    R(z)(6.312e-05 z^3 - 6.361e-05 z^2 - 6.212e-05 z + 6.262e-05)
% => x(k+4) - 3.976*x(k+3) + 5.929*x(k+2) - 3.93*x(k+1) + 0.9766*x(k) =
%    6.312e-5*r(k+3) - 6.361e-5*r(k+2) - 6.212e-5*r(k+1) + 6.262e-5*r(k)
% => x(k+4) = 3.976*x(k+3) - 5.929*x(k+2) + 3.93*x(k+1) - 0.9766*x(k) +
%    6.312e-5*r(k+3) - 6.361e-5*r(k+2) - 6.212e-5*r(k+1) + 6.262e-5*r(k)
% => x(k) = 3.976*x(k-1) - 5.929*x(k-2) + 3.93*x(k-3) - 0.9766*x(k-4) +
%    6.312e-5*r(k-1) - 6.361e-5*r(k-2) - 6.212e-5*r(k-3) + 6.262e-5*r(k-4)
%
%    H_phi(z) = P(z)/R(z) where P(z) and p(k) are for Phi
% => p(k) = 3.976*p(k-1) - 5.929*p(k-2) + 3.93*p(k-3) - 0.9766*p(k-4) +
%    -0.0001674*r(k-1) + 0.0001687*r(k-2) + 0.0001648*r(k-3) - 0.0001661*r(k-4)
%
%
% So, we get:
%  x(k) = 3.976*x(k-1) - 5.929*x(k-2) + 3.93*x(k-3) - 0.9766*x(k-4) +
%    6.312e-5*r(k-1) - 6.361e-5*r(k-2) - 6.212e-5*r(k-3) + 6.262e-5*r(k-4)
%  p(k) = 3.976*p(k-1) - 5.929*p(k-2) + 3.93*p(k-3) - 0.9766*p(k-4) +
%    -0.0001674*r(k-1) + 0.0001687*r(k-2) + 0.0001648*r(k-3) - 0.0001661*r(k-4)
%

%
% Using method from pages 49-50 of the textbook:
% http://www.me.unm.edu/~starr/teaching/me581/textbook.pdf
%
%u=zeros(3);
%y=zeros(3);
%u(0) = adc_in(); % Read input u(k).
%y(0) = a1*y(1) + a2*y(2) + b0*u(0) + b1*u(1) + b2*u(2); % compute y(k).
%dac_out(y(0)); % Write output y(k) to D/A converter.

%k = zeros(5);
%p = zeros(5);
%r = 0.0*ones(5); % Command cart to the middle, later use knob input

% Read initial values

% TODO I think this is how I come up with the estimates of the states, but then
% I also need to be doing something like this to come up with the new control
% outputs. I just need to look at the two discrete state-spaces and write out the
% difference equations from those.

% Compute new estimates
%x(1) =   3.976*x(2) -     5.929*x(3) +      3.93*x(4) -    0.9766*x(5) + \
%	  6.312e-5*r(2) -  6.361e-5*r(3) -  6.212e-5*r(4) +  6.262e-5*r(5);
%p(1) =   3.976*p(2) -     5.929*p(3) +      3.93*p(4) -    0.9766*p(5) + \
%	-0.0001674*r(2) + 0.0001687*r(3) + 0.0001648*r(4) - 0.0001661*r(5);

% Propagate variables backwards for next sample
%x(5) = x(4); x(4) = x(3); x(3) = x(2); x(2) = x(1);
%p(5) = p(4); p(4) = p(3); p(3) = p(2); p(2) = p(1);

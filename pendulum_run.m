%
% Control the inverted pendulum
%

close all;

% load values from pendulum.m
pendulum;

% Load Octave packages we'll be using
% Note: must do this below loading pendulum.m, since it clears everything
pkg load sockets control;
% load ctrlbox comm functions
ctrlbox;

% Let's run this controller.
disp('');
disp('');

try
    % Reset the board, then try connecting
    disp('Set the cart in the middle, pendulum down and stationary.')
    disp('Then click the "PROG" reset button on the FPGA board')
    disp('Finally, rotate pendulum CCW to vertical.')
    disp('');
    input('Press Enter to continue.');

    % Connection must first be established with the ctrlbox.
    ctrlbox_init();
    disp('Connected to ctrlbox');

    % Send sample period
    period = 1000000./f;
    ctrlbox_send(0,0,period);
    disp('Finished sending sample rate');

    % Initial estimates are zero
    eststate = zeros(4, 2);
    u = zeros(3);
    y = zeros(3);
    input = 0;
    estoutput = zeros(2,1);
    control_output = 0;

    % Plot output
    estoutputhistory = zeros(1,4);

    % Only run for a certain time
    cnt = 0;
    maxTime = 300; % seconds
    maxCnt = maxTime*f;

    % Let's get rid of the first bit of data so that it doesn't instantly error
    % saying the angle is too great
    for n = 1:10
        rdata = ctrlbox_recv();
    end

    %while true
    while cnt < maxCnt
        ++cnt;
        % Read encoder values
        %
        % [long_pend_angle, short_pend_angle, motor_shaft_angle, knob_angle]
        % Pendulum and motor shaft angles are 4096 counts/rev
        rdata = ctrlbox_recv();
        long_pend_angle = rdata(1)*2*pi/4096+pi;
        short_pend_angle = rdata(2)*2*pi/4096+pi;
        motor_shaft_angle = rdata(3)*2*pi/4096;
        motor_position = motor_shaft_angle*rd;
        knob_angle = rdata(4);

        % Checks for safety
        if cnt == 0 && (
            abs(motor_position) > pstart || ...
            abs(long_pend_angle) > astart)
           printf('Exiting. Invalid starting position or angle. Motor position: %f, Pendulum angle: %f\n',
               motor_position, long_pend_angle);
           break
        end

        if abs(motor_position) > Maxpos || ...
           abs(long_pend_angle) > Maxangle
           printf('Exiting. Invalid position or angle. Motor position: %f, Pendulum angle: %f\n',
               motor_position, long_pend_angle);
           break
        end

        % Set input
        %r = rdata(4)*50;
        r = 0.0;

        % Estimation
        y = [motor_position; long_pend_angle];
        uvec = [control_output; y-estoutput];
        eststate(:,1) = sys_est_only_d.a*eststate(:,2) + sys_est_only_d.b*uvec;
        estoutput = sys_est_only_d.c*eststate(:,2) + sys_est_only_d.d*uvec;
        eststate(:,2) = eststate(:,1);
        control_output = r*Nbar - K*eststate(:,1);

        % Test, use measurements in control even though we have a full-order observer
        %control_output = r*Nbar - K*[y(1); eststate(2,1); y(2); eststate(4,1)];

        % From book:
        %control_output = r*Nbar - K*(eststate(:,1) - [1;0;0;0]*r);

        % Write pwm values and enable motor
        %
        % Properly put between the +/- bounds of the PWM output.
        %
        % Scale properly, or at least in a way that appears to work.
        % Experimentally determined 2000 to work, which is close to 32768/18
        % (or 32768/20?) which came from past years' code, probably because 18
        % or 20 is the max voltage, so we're converting so the DAC will output
        % the desired voltage.
        pwm = min(max(-control_output*32768/Maxvoltage,-32767),32767);
        ctrlbox_send(pwm, 1, 0);

        % Force matlab to check for interrupts and flush event queue
        drawnow;
        fflush(stdout);

        % Save data
        printf('Long: %f, Short: %f, Motor: %f, Knob: %f, PWM: %f\n',
            long_pend_angle, short_pend_angle, motor_shaft_angle, knob_angle, pwm);
        estoutputhistory = [estoutputhistory;[estoutput' motor_position long_pend_angle]];
    end

    % Send zero
    pwm = 0;
    ctrlbox_send(pwm, 1, 0);

    drawnow;
catch
    % If something failed, display error
    disp('Exiting...');
    disp(lasterror.message);
end

% Plot estimated states
t = 0:T:(size(estoutputhistory,1)-1)/f;
figure;
[AX,H1,H2] = plotyy(t,estoutputhistory(:,1),t,estoutputhistory(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','cart position (m)');
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)');
title('State estimates');
figure;
[AX,H1,H2] = plotyy(t,estoutputhistory(:,3),t,estoutputhistory(:,4),'plot');
set(get(AX(1),'Ylabel'),'String','measured cart position (m)');
set(get(AX(2),'Ylabel'),'String','measured pendulum angle (radians)');
title('Measured values');

% disable motor and disconnect
ctrlbox_shutdown();

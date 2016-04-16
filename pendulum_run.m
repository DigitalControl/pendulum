%
% Control the inverted pendulum
%

% Load Octave packages we'll be using
pkg load sockets control;
% load ctrlbox comm functions
ctrlbox;

%
% TODO load from pendulum.m
%

% Let's run this controller.
disp('');
disp('');

try
    % Reset the board, then try connecting
    disp('Set the cart in the middle, pendulum down and stationary.')
    disp('Then click the "PROG" reset button on the FPGA board')
    disp('');
    input('Press Enter to continue.');

    % Connection must first be established with the ctrlbox.
    ctrlbox_init();
    disp('Connected to ctrlbox');

    % Send sample period
    period = 1000000./f;
    ctrlbox_send(0,0,period);
    disp('Finished sending sample rate');

    while true
        % Read encoder values
        %
        % [long_pend_angle, short_pend_angle, motor_shaft_angle, knob_angle]
        % Pendulum and motor shaft angles are 4096 counts/rev
        rdata = ctrlbox_recv();
        long_pend_angle = rdata(1);
        short_pend_angle = rdata(2);
        motor_shaft_angle = rdata(3);
        knob_angle = rdata(4);

        % Pwm values are +/-32767 for full scale motor voltage.
        % A positive pwm value causes the carriage to move in +X direction.
        % Encoder data is 32 bit raw counts, clockwise is positive.
        % Encoders are zeroed when "PROG" button is pressed on Spartan3 board.
        % Leave long pendulum hanging down when zeroing, then rotate up 180 degrees
        % when starting.  Subtract pi/2 from long_pend_angle to remove offset.
        %
        % Example control law pwm generation
        %pwm = (30000 * sin(6.28*8*double(c)/min(cnt,1000))) - rdata(1);
        pwm = rdata(4)*50;

        % Write pwm values and enable motor
        ctrlbox_send(pwm, 1, 0);

        % Force matlab to check for interrupts and flush event queue
        drawnow;
        fflush(stdout);

        % Save data
        printf('Long: %f, Short: %f, Motor: %f, Knob: %f, PWM: %f\n',
            long_pend_angle, short_pend_angle, motor_shaft_angle, knob_angle, pwm);
    end
    drawnow;
catch
    % If something failed, display error
    disp('Exiting...');
    disp(lasterror.message);
end

% disable motor and disconnect
ctrlbox_shutdown();

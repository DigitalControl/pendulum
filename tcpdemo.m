% tcpdemo4 - demo script to communicate with new pendulum interface
%
% Syntax:
%    tcpdemo4
%
% cnt is number of packets to receive
% srate is sample rate in Hz
%
% communication with the new ctrlbox is by the ctrlbox.m library
%
% receive data:
%    [long_pend_angle, short_pend_angle, motor_shaft_angle, knob_angle]
%    pendulum and motor shaft angles are 4096 counts/rev
%
% Pwm values are +/-32767 for full scale motor voltage.
% A positive pwm value causes the carriage to move in +X direction.
% Encoder data is 32 bit raw counts, clockwise is positive.
% Encoders are zeroed when "PROG" button is pressed on Spartan3 board.
% Leave long pendulum hanging down when zeroing, then rotate up 180 degrees
% when starting.  Subtract pi/2 from long_pend_angle to remove offset.
%
% Protocol is send pwm value, then read encoders.  Sampling in ctrlbox is
% done according to a hardware sample clock.
%
pkg load sockets control;
ctrlbox;            % load ctrlbox comm functions

cnt=4000;            % number of times through loop
srate = 400;            % sample rate in Hz

store = zeros(cnt,5);
rdata = [0,0,0,0];        % receive data

% Connection must first be established with the ctrlbox.

try
    ctrlbox_init();

    disp('finished init');

    % send sample period
    period = 1000000./srate;

    ctrlbox_send(0,0,period);

    disp('finished send');

    x = 1:cnt;
    tic;
        for c=1:cnt
        % read encoder values
        rdata = ctrlbox_recv();

        %
        % your control law goes here
        %

        % example control law pwm generation
        pwm = (30000 * sin(6.28*8*double(c)/min(cnt,1000))) - rdata(1);

        % write pwm values and enable motor
        ctrlbox_send(pwm, 1, 0);

        % force matlab to check for interrupts and flush event queue
        drawnow;

        % save data
        store(c,:) = [rdata,pwm];
    end
    runtime = toc;
    fprintf('transactions=%d seconds=%d transactions/sec=%f\n',
        c, runtime, c/runtime);
    drawnow;

catch
    % if something failed, display error and loop count
    fprintf('c=%d\n',c);
    disp(lasterror.message);
end
% disable motor and disconnect
ctrlbox_shutdown();
plot(store)

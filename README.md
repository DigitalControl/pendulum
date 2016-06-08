# Digital Control of Inverted Pendulum
### Connecting to the Pendulum
This is using Walla Walla University's inverted pendulum up in the digital or
analog labs in Kretschmar.

#### Connect to board
Plug into ethernet, set IPv4 settings:

    IP: 169.254.0.2
    Netmask: 255.255.255.0
    Gateway: 169.254.0.2

Set the ethernet MAC address to be (Identity --> Cloned Address):

    00:13:F7:1A:79:FA

Check "Use this connection only for resources on its network"

Reset the board. Then connect to this with NetworkManager if you're using it.

#### Install Octave libraries
On Arch Linux, you can install from the AUR, e.g.:

    spinach -i octave-control octave-sockets

or, in Octave, download the [control](http://octave.sourceforge.net/control/)
and [sockets](http://octave.sourceforge.net/sockets/) packages and then
install them:

    pkg install control-3.0.0.tar.gz sockets-1.2.0.tar.gz

#### Testing Pendulum Connnection
There's a _tcpdemo_ script that'll make the motor move back and forth and plot
the data it received back over the network connection.

### Pendulum Control
Run in simulation:

    octave
    >> pendulum            # 4th order, 30 V, correct rd
    >> pendulum_6th        # 6th order, 30 V, correct rd
    >> pendulum_wrong      # 4th order, 20 V, wrong rd (rd^2)
    >> pendulum_spring     # 8th order, not observable or controllable

In these non-_\_run_ files, at the top is a _plotAll_ variable that if set to true will generate all the plots. In the 6th order one, there's also a _bothPendulums_ variable which if set to true will balance both pendulums rather than the long one up and the small one down.

Run on the pendulum:

    octave
    >> pendulum_run        # 4th order, 30 V, correct rd
    >> pendulum_6th_run    # 6th order, 30 V, correct rd
    >> pendulum_wrong_run  # 4th order, 20 V, wrong rd (rd^2)

In the _\_run_ files, there's a _saveData_ variable which if set to true will save data while the system is running and then stop balancing the pendulum after a certain number of seconds (set by _maxTime_) and generate plots of the measurements and estimates.

### Report
The report documenting the derivations and results is [Pendulum\_project\_report.pdf](http://github.com/DigitalControl/pendulum/blob/master/Pendulum_project_report.pdf).

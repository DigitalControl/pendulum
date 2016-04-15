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

or, in Octave:

    pkg install control sockets

#### Testing Pendulum Connnection
There's a _tcpdemo_ script that'll make the motor move back and forth and plot
the data it received back over the network connection.

### Pendulum Control
_yet to be implemented_

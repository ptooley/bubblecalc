# Bubblecalc
BubbleCalc is used to calculate the motion of electrons in an idealised wakefield
bubble. Electron bunch and bubble parameters are specified via a configuration
file. Output is in hdf5 format with SI units of Orisis units.

# Usage
BubbleCalc has no real command line interface, but accepts a
single .ini formatted input file which contains all options. This is a deliberate
design decision to encourage documentation and reproducibility of research.
```
user@machine $ bubblecalc −h
Usage: bubblecalc [−h] <configfile>
Allowed Options:
−h [−−help]  Print this help message
<configfile> Configuration .ini file
```

# Configuration File
The configuration file format is .ini, with key = value pairs. The full set of available options are shown below, note that not all fields are required for a valid config, default options are shown where applicable:

```ini
; Not all fields are required for a valid config, default options are shown where applicable

; Output section contains file handling options
[Output]
; number of timesteps to split trajectory into
nsteps = [no default, integer]

; output filename
filename = [no default, string - must be valid filename]

; make output units compatible with radt (Author use)
radt_mode = false

; Bubble section contains global physics related parameters
[Bubble]
; laser wavelength providing reference for critical density (units of meters)
lambda0 = [no default, real number]

; plasma density as a fraction of the critical density for given laser wavelength
eta2 = [no default, real number]

; Bunch sections specify properties of individual bunches
; A single file can contain an arbitrary number of bunch sections
; all bunch sections should have a title [Bunch*], e.g
[Bunch1]
; All units are SI, laser propagation is in the positive x direction
; Input coordinate system is referenced relative to the centre of the bubble
; this means typically bunches will be specified with negative x position

; number of particles in the bunch
npart = [no default, integer]

; total bunch charge, if not set assume each particle represents one electron
q = [no default, integer]

; statistical distributions are set as follows, example of x position:
; "normal" normal distribution with centre x, standard deviation dx
; "sin2" sin-squared distribution with centre x and width dx
; "constant" all values x, dx ignored
; "linspaced" uniformly spaced over the range x-dx/2 to x+dx/2
; "" - if unset and dx given, normal distribution, otherwise constant
;xdist = ""

; average x position, defaults to 0
;x = 0.0

; width of x distribution, defaults to 0
;dx = 0.0

; These options are similarly set for the y, z, positions, px, py, pz momentum distributions
; and start time t:

;ydist = ""
;y = 0.0
;dy = 0.0

;zdist = ""
;z = 0.0
;dz = 0.0

;pxdist = ""
;px = 0.0
;dpx = 0.0

;pydist = ""
;py = 0.0
;dpy = 0.0

;pzdist = ""
;pz = 0.0
;dpz = 0.0

;tdist = ""
;t = 0.0
;dt = 0.0
```
# Output File Structure

The output HDF5 file will contain one group per electron specified in the input
file. Each electron group contains time series datasets for the electron motion in
either SI units or Osiris units as requested by the user.

#Installation
Source code is available from the github repository: https://github.com/ptooley/bubblecalc.git

Prerequisites for installation are:
• A C++11 capable compiler (e.g gcc 6 or greater)
• Cmake version 3.6 or later
• Boost version 1.58 or later
• HDF5 (version 1.10 or later recommended)
Installation is then performed by checking out the code, running cmake and
then make:

```
user@machine $ git clone https://github.com/ptooley/bubblecalc.git
user@machine $ cd bubblecalc
user@machine $ cmake .
user@machine $ make
```
BubbleCalc can then be found in the ./bin/ directory.

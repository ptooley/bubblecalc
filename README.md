# Bubblecalc


# Configuration File

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

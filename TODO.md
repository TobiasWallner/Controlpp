
TODO:
=====

- bode shall use frequencies in rad/s instead of hz.
- make type safe types for:
    - rad/s
    - hz
    - absolute magnitudes
    - magnitudes in dB

- Method to interpolate timestep data

- make a `discretize_zoh()` and a `discretize_tustin()` function for discretisation.

Hybrid & Dynamic Storage
--------------
Max storage on stack with dynamic size
- Adapt Polynomial to allow hybrid vectors
- Adapt transfer functions to use hybrid vectors
- Adapt state space to allow hybrid matrices
- Adapt math functions to allow for hybrid vectors/matrices

State Space Arithmetic
----------------------
- Add
- Subtract
- Multiply
- Divide

Evaluate Transfer Functions
---------------------------
- Bodesweep
- Stepresponse
- Complex evaluation
- Zeros and Poles --> find roots of polynomials

Controllers
-----------
- LQR (Linear Quadratic Regulator)
- LQG (LQR + Kalman)
- H_2 (in progress) <----------------------------------
- H_inf

Filters
-------
- Saturation
- SoftSaturation

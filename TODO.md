
TODO:
=====

- bode shall use frequencies in rad/s instead of hz.
    - offer convenient conversion functions for scalars and vectors:
        - radps_to_hz
        - hz_to_radps
        - abs_to_dB
        - dB_to_abs

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

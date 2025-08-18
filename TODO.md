
TODO:
- use 'z^{-1}' instead of 'z' for discrete transfer functions.

- larger tests of controllers designed with `tf::s`

Naming Style
------------
- All functions and Methods in `snake_case`
- Name controllers like: `PController`, `PControl` or `controller::P` instead of just `P`
  >Looking at this I think I prefer `PControl`

Hybrid Storage
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

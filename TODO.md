
TODO:
- Check what is going wrong in the Clang Build

Naming Style
------------
- Functions and Methods in `snake_case`

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

Filters
-------
- Saturation
- SoftSaturation

Time Variable Controllers:
--------------------------
- P
- I 
  - Saturated
  - RateLimited
- D
  - Tamed O1
  - Tamed O2
- PT1
- PT2
- PID
  - Saturated I
  - RateLimited I
  - Tamed D O1
  - Tamed D O2
- PTID
  - Saturated I
  - RateLimited I
  - Tamed D O1
  - Tamed D O2
- LeadLag
- Notch
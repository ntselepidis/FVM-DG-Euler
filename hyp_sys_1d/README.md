# FVM and DG for the Euler equations in 1D
The exercise requests implementing a large number of different schemes, which
when combined lead to a very large number of different solvers. Doing this in
finite time and in a somewhat maintainable manner requires somewhat complex
code. This document will help you navigate the code.

The idea is quite simple: split up the FVM/DG into its individual components and
implement each of these and then assemble the complete method by aggregation.

The different components for FVM are:
  
  * model, you have to implement the Euler equations, Burgers equation is given as an example.

  * reconstruction
    - piecewise constant
    - piecewise linear

  * slope-limiters

  * numerical fluxes

  * boundary conditions

  * one-step time integration for the semi-discrete version of
      du/dt = - div(f(u)) =: L(u)

  * I/O, input is done for you. The main structure is provided for output, you need to implement which quantity of the       system you want to write to file.

Other than the above-mentioned components, excluding slope-limiters, DG additionally requires the following:

  * polynomial basis.
  
  * A handler, which provides different routines specific to DG that would be used.
  
  * A limiting procedure which is applied to the solution coefficients.

FVM components are combined as follows:

  * flux-loop which computes du/dt (or L(u) = div(u) in the notation above). This
    relies on the reconstruction (including the slope-limiter) and numerical flux.

  * the one-step time integration uses the boundary condition and flux-loop.

  * the loop which chains multiple steps relies on the single step update and the I/O.

DG components are combined as follows:

  * Update-loop which computes du/dt (or L(u) = div(u) in the notation above). This
    relies on the numerical flux, polynomial basis and DG handler.

  * the one-step time integration uses the update-loop, boundary condition and DG limiting.

  * the loop which chains multiple steps relies on the single step update and the I/O.

## Testing
It is very beneficial to write short tests to ensure your code does what you
think it's doing. In a first step this does not need to be rigorous even simple
tests are beneficial.

You can use writing tests as a means of debugging, i.e. if you want to know if a
piece of code does X, you write a test to check that it does X.

## How to get access to the simulation time?
If you at any point need access to the current time in the simulation or the
current time-step, e.g. Lax-Friedrichs flux. You will need to store a shared
pointer to the `SimulationTime` object used by the `TimeLoop`. The
`JSONSnapshotWriter` is an example of this technique.

## Where to start?
First look at `TimeLoop`. Look at what other (polymorphic) objects it uses. Then
look at each of these to get a sense of the overall architecture. Maybe you want
to sketch it out.

Now you can start the exercise, i.e. implement a numerical flux. Note that
solving the exercise in the intended order is beneficial. However, you do not
need to implement the whole sub-exercise at once, you can choose to only do one
numerical flux and do the others later.

## Why do we need to 'register' our components?
The idea is to be able to choose the numerical method at run-time through a
configuration file called `config.json`. This requires that each component has
a name e.g. `"ssp2"` which can be used to create an object of type `SSP2` (as
opposed to `ForwardEuler`).

This mapping of string identifier to class we implement in functions `make_*`.
We choose the simplest version of implementing these namely a bunch of
if-statements, e.g.
```
    if(key == "ForwardEuler") {
        return std::make_shared<ForwardEuler>(...);
    }

    if(key == "ssp2") {
        return std::make_shared<SSP2>(...);
    }

    throw std::runtime_error("...");
```

Finally to reduce the repetition and clutter in the code above, we provide macros
which would turn the above into:
```
    REGISTER_RUNGE_KUTTA("forward_euler", ForwardEuler);
    REGISTER_RUNGE_KUTTA("ssp2", SSP2);

    throw std::runtime_error("...");
```

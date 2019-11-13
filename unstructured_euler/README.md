# NACA airfoil simulation
In this exercise you'll simulate the steady flow over an airfoil!

## Coding rules
This is not an exercise in creating versatile code. You only need to implement
one algorithm, and you may assume that you will never be asked for a second
version of the algorithm. This allows you to save a lot of boiler plate code,
and allows you to make simplifying assumptions.

This does not mean that you do not need to write comprehensible code. It only
means you do not have to write particularly general code.

Remember, the goal is to simulate an airfoil with a FVM with piecewise linear
reconstruction of the primitive variables and HLLC numerical flux.

### Order of implementation
The easiest way to implement this algorithm is to start with a first order
solver. It requires:

  - the HLLC (or Rusanov) numerical flux
  - outflow flux boundary conditions
  - CFL condition

You can then run the vortex simulation an check convergence.

Step two, implement reconstruction in conserved variables. The only change
required is computing the trace values. You again check convergence rates for
the vortex problem.

Step three, implement reconstruction in primitive variables. This is a simple
modification of your existing code; and again, you check convergence rates for
the vortex problem.

### Available libraries
You have Eigen, JSON and libfmt at your disposal. If you want to use GTest, it's
also been set up, you need to add `-DHAS_GTEST=1` to your `cmake` command.
However, take care to implement all tests in one compilation unit, since
otherwise you'll get many link time errors about duplicate definitions.

## Compilation
The performance of the code is *highly* reliant on high levels of optimization.
For this purpose we've added an additional build type called `FastDebug` it is
`O2 -g`, i.e. you get high levels of optimization and debug symbols,
additionally it does not set `NDEBUG=1` (which deactivates any `assert` in your
code and Eigen, i.e. you wont have any bounds checks anymore). I would recommend
you use `FastDebug` unless you can't figure out what's happening in the debugger
(use `Debug`) or once you've completely debugged you code and simply want to run
the simulation (use `Release`).

Note:
  - `FastDebug` is 100x faster than `Debug` on my laptop.
  - `Release` is 2x faster than `FastDebug` on my laptop.

Therefore, you will need at least two build types, `Debug` and
`FastDebug`/`Release`. If you are using an editor and a shell you simply
generate two build folders:

    mkdir build-debug && cd build-debug && cmake -DCMAKE_BUILD_TYPE=Debug ..
    cd ..
    mkdir build-fast-debug && cd build-fast-debug && cmake -DCMAKE_BUILD_TYPE=FastDebug ..

If you're using an IDE you can set up multiple build types. (In CLion you find
this under Settings -> Build, Execution, Deployment -> CMake; remember to add
the build type as a flag to CMake in the box "CMake options")

### OpenMP
You can enable OpenMP by using passing `-DHAS_OPENMP=1` to `cmake`. This will
enable OpenMP support in the compiler. You can then parallelize the loops over
all cells of the grid using the required OpenMP pragma. This can buy you more
speed up :)

### Profiling
If you feel you need to profile you code, we've set up support of gperftools
[1]. This is a sampling profiler you can use it by adding `-DHAS_GPERF=1` to
your `cmake` command. You need to relink your binary and execute it as follows:

    CPUPROFILE=vortex.prof ./vortex

This will generate a profile `vortex.prof` which you can analyse using `pprof`,
e.g.

    pprof ./naca_airfoil naca.prof
    # or accumulated version
    pprof -text -cum ./naca_airfoil naca.prof


[1]: https://github.com/gperftools/gperftools

### Static code analysis
To enable `clang-tidy` you need to add `-DHAS_CLANG_TIDY=1` to your `cmake`
command.

## Numerical Experiments
There are two numerical experiments ready for you.

### Isentropic Vortex
This initial condition is coded up in `vortex.cpp`. It is a smooth test
case of a flow spiraling out from the center of the domain. We've included a
reference solution with the templates, you can find it in `data/vortex`.

You can use

    python ../compute_convergence.py vortex_{3..5} ../data/vortex/vortex_7

to compute the error and convergence rate (replace `{3..5}` as needed, but do
not skip any refinements level in the middle).

You expect to see some convergence rates of 0.8 - 1.0 for your first order
solver; and 1.7 - 2.0 for you second order solver. Note, that the small grids
are all preasymptotic, therefore they do not yet converge at the expected rate,
i.e. they can sub- or super-convergent. The final grid is a bit too close to the
reference solution. So you can again accept odd behaviour.

Note, that you will see something call grid imprinting. That is, you might see
features of the grid clearly visible in the numerical solution.

You can plot the grid and solution with

    python ../plot_on_mesh.py vortex_4


### NACA airfoil
The second numerical experiments is the simulation of an airfoil. It is coded up
in `naca_airfoil.cpp`. You should only (seriously) start running this test after
you've measured the expected convergence rates for the vortex.

You can modify the angle of attack and Mach number of the initial conditions and
plot the results using `plot_on_mesh.py`.

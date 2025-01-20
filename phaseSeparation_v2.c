/**
# Coupled reaction--diffusion equations

The [Brusselator](http://en.wikipedia.org/wiki/Brusselator) is a
theoretical model for a type of autocatalytic reaction. The
Brusselator model was proposed by Ilya Prigogine and his collaborators
at the Free University of Brussels.

Two chemical compounds with concentrations $C_1$ and $C_2$ interact
according to the coupled reaction--diffusion equations:
$$
\partial_t C_1 = \nabla^2 C_1 + k(ka - (kb + 1)C_1 + C_1^2 C_2)
$$
$$
\partial_t C_2 = D \nabla^2 C_2  + k(kb C_1 - C_1^2 C_2)
$$

We will use a Cartesian (muDifflti)grid, the generic time loop and the
time-implicit diffusion solver. */

#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "tension.h"
#include "diffusion.h"
/**
We need scalar fields for the concentrations. */

scalar C1[];

/**
We use the same parameters as [Pena and Perez-Garcia,
2001](/src/references.bib#pena2001) */

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

/**
## Parameters

We change the size of the domain `L0` and set the tolerance of the
implicit diffusion solver. */

u.n[left] = neumann(0.);
u.n[right] = neumann(0.);
p[left] = dirichlet(0.);
p[right] = dirichlet(0.);

int main()
{
  init_grid (256);
  L0 = 4e0;
  X0 = -L0/2.;
  TOLERANCE = 1e-4;
  rho1 = 1e0; rho2 = 1e-1;
  mu1 = 1e-1; mu2 = 1e-1;
  char comm[160];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  /**
  Here $\muDiff$ is the control parameter.  For $\muDiff > 0$ the system is
  supercritical (Hopf bifurcation). We test several values of $\muDiff$. */

  f.sigma = 1e0;
  f.tracers = {C1};
  
  run();
}

/**
## Initial conditions */

event init (i = 0)
{
  fraction(f, sq(1.) - sq(x) - sq(y));

  foreach() {
    C1[] = C1[] = f[]*(x>0. ? 1. : -1.)*noise(); 
  }
  // foreach(){
  //   u.x[] = C1[];
  // }
}

// event acceleration (i++)
// {
//   face vector av = a;
//   foreach_face(x){
//     double ff = (f[] + f[-1])/2.;
//     av.x[] += 1e1*ff*(C1[-1] + C1[])/2.*exp(-t);
//   }
//   foreach_face(y){
//     double ff = (f[] - f[0,-1])/Delta;
//     av.y[] += (fabs(x)<2e-1?1.:0.)*2e1*ff;
//   }
// }

event tracer_diffusion (i++){
  // foreach(){
  //   u.x[] = C1[];
  // }
  face vector diffusion_coef[];
  foreach_face(){
      double ff = (f[] + f[-1])/2.;
      diffusion_coef.x[] = 0.05*fm.x[]*ff;
  }
  diffusion (C1, dt, D = diffusion_coef);

  foreach(){
    C1[] *= f[];
  }

}

event output(t = 0, t += 1e-3, t <= 0.1)
{
  char dumpFile[160];
  sprintf (dumpFile, "intermediate/snapshot-%5.4f", t);
  dump (file = dumpFile);
}

event adapt(i++){    
  adapt_wavelet({f, u.x, u.y, C1}, (double[]){1e-2, 1e-1, 1e-1, 1e-1}, 11);
}
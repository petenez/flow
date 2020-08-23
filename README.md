# flow
Computational fluid dynamics with phase field crystal etc.

## float
<img src="https://github.com/petenez/flow/blob/master/float/demo/flow-10000.png" width="320"/> <img src="https://github.com/petenez/flow/blob/master/float/demo/phi-10000.png" width="320"/>

This project combines phase field crystal (PFC) with computational fluid dynamics (CFD) to model something like ice crystals growing and floating in water. See an animated simulation (https://www.youtube.com/watch?v=QUsN1aDdjC4) or try it yourself by running the scripts simulate.sh and visualize.sh in flow/float/demo. See my thesis (http://urn.fi/URN:ISBN:978-952-60-8608-8) and my GitHub repository for PFC tools (https://github.com/petenez/pfc) for further details on PFC. The CFD solver is based on the paper Real-Time Fluid Dynamics for Games by Jos Stam (available: https://pdfs.semanticscholar.org/847f/819a4ea14bd789aca8bc88e85e906cfc657c.pdf).

A partially crystalline system is modeled using a basic PFC model. The PFC density field describing the crystals is advected along the CFD flow field. The flow is also coupled to the smoothed density, where the lattice structure has been smeared out, in two ways. First, an upward force proportional to the smoothed density is applied to the flow, causing the crystals to float, the faster the larger they are. Second, the viscosity of the underlying fluid is also proportional to the smoothed density to make the flow interact with the crystals more like rigid bodies. The net flow rate is fixed to zero, resulting in equal flow up and down, and left and right.

The actual simulation code is written in C and is parallelized with OpenMP. The PFC dynamics are solved exploiting fast Fourier transforms that are computed with FFTW. The visualization tools require Java and MPI.

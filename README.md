# Car Exhaust

In car exhaust tips, is it possible to control the 
frequencies of sound leaving?  This is tested by adding fins into a blocked
side of the exhaust tips.  

## Inner Workings

This project examines the problem using a 2 dimensional approximation to 
the exhaust tips chambers.  If $\Omega$ is the domain of the chambers,
let $\Gamma_{in}$ be the input side, $\Gamma_{out}$ be the output side
and boundary by $\partial \Omega$.  Let $u$ be the pressure,
$t$ be the time, and $x \in \mathbb R^2$ be the location.  

We add a collection of fins, occupying space 
$\Omega_{fin}$.  These act as wave reflectors and break out the sound.  
By changing the spacing and location we can change the accoustic 
properties of the exhaust tips.

The problem can be summarized as follows:

The governing PDE is
$$
\frac{\partial^2 u}{\partial t^2} = \nabla^2 u, \quad 
x \in \Omega \setminus \Omega_{fin}
$$

The inbound boundary term is
$$
\left. u(x, t) \right|_{x \in \Gamma_{in}} = f(t)
$$

The outbound boundary condition is the one way wave equation
$$
\left(\frac{\partial u}{\partial t} + \frac{\partial u}{\partial n} 
\right)_{x \in \Gamma_{out}}= 0
$$
Where $n$ is the unit outward normal.

Finally everywhere else on the boundary we have a Neumann boundary
condition
$$
\frac{\partial u}{\partial n} = 0, \quad 
x \in \partial \Omega \setminus (\Gamma_{in} \cup \Gamma_{out}) \cup
\partial \Omega_{fin}
$$

The problem is solved in FEniCS using the variational formulation.  
The time stepping is done using the Newmark-$\beta$ method.

The script PDE.py is the main code for the problem.  Supporting 
codes are found in TubularTuber.py, which builds the domain and
WaveBreaker.py, which constructs the fins.

## Dependencies

The project relies on the FEniCS project libraries as well as the 
mshr scripts to mesh the domain.  It also relies on scipy, numpy,
and matplotlib.  These can be installed using the Anaconda software 
by running

\$ conda install -c conda-forge fenics mshr numpy scipy matplotlib

It is recommended to install FEnICs to a separate conda environment. 
Directions for this are located on the FEniCS project website:
https://fenicsproject.org/

## Additional Software

The data files are saved in the XDMF format.  This has an hdf5 file
as the data set and an XDMF file as the overview.  It is easiest 
to view the solution using Paraview.  Installation instructions 
can be found here: https://www.paraview.org/.

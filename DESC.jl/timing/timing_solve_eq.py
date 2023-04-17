import time 
start = time.time()
import numpy as np
from desc.equilibrium import Equilibrium
from desc.geometry import FourierRZToroidalSurface
from desc.profiles import PowerSeriesProfile
from desc.plotting import plot_1d, plot_section, plot_surfaces


surface = FourierRZToroidalSurface(
    R_lmn=[10, 1],
    modes_R=[[0, 0], [1, 0]],  # modes given as [m,n] for each coefficient
    Z_lmn=[0, -1],
    modes_Z=[[0, 0], [-1, 0]],
)

pressure = PowerSeriesProfile(params=[0, 0], modes=[0, 2])
iota = PowerSeriesProfile(params=[1, 1.5], modes=[0, 2])

eq = Equilibrium(
    surface=surface,
    pressure=pressure,
    iota=iota,
    Psi=1.0,  # flux (in Webers) within the last closed flux surface
    NFP=1,  # number of field periods
    L=7,  # radial spectral resolution
    M=7,  # poloidal spectral resolution
    N=0,  # toroidal spectral resolution (axisymmetric case, so we don't need any toroidal modes)
    L_grid=12,  # real space radial resolution, slightly oversampled
    M_grid=12,  # real space poloidal resolution, slightly oversampled
    N_grid=0,  # real space toroidal resolution (axisymmetric, so we don't need any grid points toroidally)
    sym=True,  # explicitly enforce stellarator symmetry
)

from desc.optimize import Optimizer

# create an Optimizer object using the lsq-exact optimizer
optimizer = Optimizer("lsq-exact")

from desc.objectives import (
    get_fixed_boundary_constraints,
    ObjectiveFunction,
    FixBoundaryR,
    FixBoundaryZ,
    FixPressure,
    FixIota,
    FixPsi,
    ForceBalance,
)

constraints = (
    FixBoundaryR(),  # enforce fixed  LCFS for R
    FixBoundaryZ(),  # enforce fixed  LCFS for R
    FixPressure(),  # enforce that the pressure profile stay fixed
    FixIota(),  # enforce that the rotational transform profile stay fixed
    FixPsi(),  # enforce that the enclosed toroidal stay fixed
)
# choose the objectives to be ForceBalance(), which is a wrapper function for RadialForceBalance() and HelicalForceBalance()
objectives = ForceBalance()
# the ObjectiveFunction object which we can pass to the eq.solve method
obj = ObjectiveFunction(objectives=objectives)

constraints2 = get_fixed_boundary_constraints()

eq.solve(
    verbose=3,
    ftol=1e-8,
    maxiter=50,
    constraints=constraints,
    optimizer=optimizer,
    objective=obj,
)


end = time.time()

print("Overall time: ", end - start)


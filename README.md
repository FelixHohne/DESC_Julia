# DESC_Julia

This is a PyCall-based Julia wrapper, with GPU acceleration, for the [DESC library](https://github.com/PlasmaControl/DESC). It supports the API as provided by the API Documentation, available [here](https://desc-docs.readthedocs.io/en/stable/api.html).  

### Nvidia GPU Installation Instructions

Note: these instructions assume CUDA 11. On Cornell’s G2 cluster, one can verify CUDA 11 is present via (`cd /usr/local/` and checking cuda-11.2 is present). 

In particular, note that DESC does not always support the latest version of JAX, and thus older versions of GPU-accelerated versions of JAX may need to be installed, which may in turn require lower versions of JaxLib, as well as CUDA and CuDNN. 

1. Run `git clone https://github.com/PlasmaControl/DESC`
2. Open the requirements.txt file in the cloned DESC folder, and remove line 3 (jax[cpu] >= 0.2.11, <= 0.4.1). DESC should now have no Jax dependencies. This is necessary as currently DESC will remove any versions of JAX with GPU-acceleration if the line jax[cpu] is present in the README. 
3. Create a new conda environment (i.e. `conda create --name desc_jl python=3.9`) (Note: for now, python = 3.11 is not supported, because at least in April 2023, some key dependencies do not yet support python 3.11.  
4. Ensure GPU access is available on the environment (i.e. for Cornell University’s G2 Cluster, run `srun -N 1 --gres=gpu:1 --mem=10000 --time=12:00:00 --pty bash`)
4. Activate the conda environment (`conda activate desc_jl`)
5. cd into the DESC/ directory 
6. Run `pip install jax==0.3.24 jaxlib==0.3.24+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`. Note that this step may fail if the current environment does not have CUDA properly configured (i.e. it is run on a login server)   
7. Verify the successful GPU installation of Jax by opening the python interpreter on the command line, then running `import jax`, followed by `jax.default_backend()`. This should report “gpu”. 
8. Run `pip install -e .`in the root directory of DESC. This will install Python DESC as a python package. 
9. cd out of the DESC/ directory. 
10. Run `git clone https://github.com/FelixHohne/DESC_Julia`. 
11. cd into DESC_Julia/DESC.jl/
11. Open Julia on command line. 
12. Run `using Pkg`, followed by `]`
13. Run `activate .` The left side of the command prompt should now be (DESC) pkg. 
14. Run `build` 
15. Run `test DESC`. In the command line output, it should report the GPU being used. 

### CPU Installation (tested on M1 Mac)
1. Create a new conda environment (i.e. `conda create --name desc_jl python=3.9`) (Note: for now, python = 3.11 is not supported, because at least in April 2023, key dependency numba does not yet support python 3.11.  
2. Activate the conda environment (`conda activate desc_jl`) 
3. Install Python DESC: `pip install desc-opt`. We can install DESC normally here as the jax[cpu] dependencies do not cause any issues. 
4. Run `git clone https://github.com/FelixHohne/DESC_Julia`. 
5. cd into DESC_Julia/DESC.jl/
6. Open Julia on command line. 
7. Run `using Pkg`, followed by `]`
8. Run `activate .` The left side of the command prompt should now be (DESC) pkg. 
9. Run `build` 
10. Run `test DESC`. It should report CPU is being used. 

### To use DESC.jl:
1. Open Julia in the terminal
2. cd into `DESC/DESC.jl`
3. Run `using Pkg`
4. Run `Pkg.activate(“.”)`
5. Run `using DESC`


### General Usage 

The library contains initialization functions in Julia for DESC objects, which use PyCall to wrap Python objects. The names of these methods are the same as DESC, as long as the naming does not result in conflicts. For example, Python DESC constructs Equilibriums at `desc.equilibrium.Equilibrium`, which is implemented here at `descEquilibrium.Equilibrium`. In a similar vein, Python DESC's `desc.plotting.plot_boozer_surface` is available at `descPlot.plot_boozer_surface`. 

Through PyCall, methods on DESC objects can be applied directly. Additional Julia functions are provided for class method functions, as well as methods which utilize numpy arrays. PyCall does not normally convert arrays to C contiguous style, leading to possible errors in the underlying DESC Python code. 

For example, Equilibrium perturb functionality is provided as a method call on a DESC Equilibrium object, i.e. eq.perturb(...). If `julia_eq` is a Equilibrium object constructed by this library, then likewise, one can call `julia_eq.perturb(...)` passing in Julia objects. However, if numpy arrays are passed into `julia_eq`, these may result in errors. As such, `descEquilibrium.equilibrium_perturb(julia_eq, ...)` provides similar functionality, where the arguments to this function are wrapped. 

### Example Program:
```
DESC.use_gpu_if_available(); # Checks if GPU is available using JAX, defaults to CPU is not available. 
# Default arguments match Python DESC; thus, only need to pass in arguments that are non-default. 
surface = DESC.FourierRZToroidalSurface(
        R_lmn=[10, 1],
        modes_R=[[0, 0], [1, 0]], 
        Z_lmn=[0, -1],
        modes_Z=[[0, 0], [-1, 0]],
);

knots = LinRange(0, 1, 20);
pressure_values = zeros(20);
iota_values = 1 .+ 1.5 * (knots.^2);

pressure = DESC.SplineProfile(
	values = pressure_values, knots = knots
);

iota = DESC.SplineProfile(
	values = iota_values, 
	knots = knots
);

# Example of how to return proper Julia objects from DESC_Julia
julia_iota = convert(Array{Float64}, iota.params);

eq = DESC.Equilibrium(
    surface=surface,
    pressure=pressure,
    iota=iota,
    Psi=1.0,  
    NFP=1,  
    L=6,  
    M=6,  
    N=0, 
    L_grid=12,  
    M_grid=9, 
    N_grid=0,  
    sym=true
);

optimizer = DESC.Optimizer("lsq-exact");

constraints = (
        DESC.FixBoundaryR(), 
        DESC.FixBoundaryZ(), 
        DESC.FixPressure(), 
        DESC.FixIota(), 
        DESC.FixPsi()
);

objectives = DESC.ForceBalance();

obj = DESC.ObjectiveFunction(objectives);

eq.solve(
    verbose=2, ftol=1e-8, objective=obj, optimizer=optimizer, constraints=constraints
);

plot = DESC.plot_section(
  eq, 
  "|F|", 
  norm_F=true, log=true,
  figsize=(12, 12), title_font_size = 25,
  xlabel_fontsize = 20, ylabel_fontsize = 20, 
);

```

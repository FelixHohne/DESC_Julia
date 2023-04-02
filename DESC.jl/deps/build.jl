using PyCall
import Pkg

ENV["PYTHON"] = ENV["CONDA_PREFIX"] * "/bin/python"
println(ENV["PYTHON"])
Pkg.build("PyCall")
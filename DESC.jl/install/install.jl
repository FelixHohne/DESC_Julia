function desc_install()
    ENV["PYTHON"] = ENV["CONDA_PREFIX"] * "/bin/python"
    println("Our conda environment is: ", ENV["PYTHON"])
    Pkg.build()
end


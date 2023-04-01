function pygranso_install()
    ENV["PYTHON"] = "/opt/miniconda3/envs/arm_ml/bin/python"  
    Pkg.build() 
    
end

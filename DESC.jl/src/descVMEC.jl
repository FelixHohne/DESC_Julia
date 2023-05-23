using PyCall

#= 
Functions implement VMEC IO API
=#


function vmecio_compute_coord_surfaces(
    equil, 
    vmec_data; 
    Nr = 10, 
    Nt = 8, 
    Nz = nothing,  
    kwargs...
    )
    
    kwargs_dict = Dict(pairs(kwargs)) 

    py"""
    import desc 
    import desc.vmec

    def vmec():
        return desc.vmec.VMECIO.compute_coord_surfaces(
           $equil, 
           $vmec_data, 
           Nr = $Nr, 
           Nt = $Nt, 
           Nz = $Nz, 
           **$kwargs_dict
        )
    """
    output = py"vmec"()
end 

function vmecio_compute_theta_coords(
    lmns, 
    xm, 
    xn, 
    s, 
    theta_star, 
    zeta; 
    si = nothing
    )
    
    py"""
    import desc 
    import desc.vmec
    if isinstance($lmns, np.ndarray):
        new_lmns = np.ascontiguousarray($lmns)
        assert new_lmns.flags['C_CONTIGUOUS']
    else:
        new_lmns = $lmns
    if isinstance($xm, np.ndarray):
        new_xm = np.ascontiguousarray($xm)
        assert new_xm.flags['C_CONTIGUOUS']
    else:
        new_xm = $xm

    if isinstance($xn, np.ndarray):
        new_xn = np.ascontiguousarray($xn)
        assert new_xn.flags['C_CONTIGUOUS']
    else:
        new_xn = $xn

    if isinstance($s, np.ndarray):
        new_s = np.ascontiguousarray($s)
        assert new_s.flags['C_CONTIGUOUS']
    else:
        new_s = $s

    if isinstance($theta_star, np.ndarray):
        new_theta_star = np.ascontiguousarray($theta_star)
        assert new_theta_star.flags['C_CONTIGUOUS']
    else:
        new_theta_star = $theta_star

    if isinstance($zeta, np.ndarray):
        new_zeta = np.ascontiguousarray($zeta)
        assert new_zeta.flags['C_CONTIGUOUS']
    else:
        new_zeta = $zeta

    if isinstance($si, np.ndarray):
        new_si = np.ascontiguousarray($si)
        assert new_si.flags['C_CONTIGUOUS']
    else:
        new_si = $si
        
    def vmec():
        return desc.vmec.VMECIO.compute_theta_coords(
           new_lmns, 
           new_xm, 
           new_s, 
           new_theta_star, 
           new_zeta, 
           si = new_si
        )
    """
    output = py"vmec"()
end 


function vmecio_load(
    path; 
    L = nothing, 
    M = nothing, 
    N = nothing, 
    spectral_indexing = "ansi", 
    profile = "iota"
) 
    py"""
    import desc 
    import desc.vmec

    def vmec():
        return desc.vmec.VMECIO.load(
           $path, 
           L = $L, 
           M = $M, 
           N = $N, 
           spectral_indexing = $spectral_indexing, 
           profile = $profile
        )
    """
    output = py"vmec"()
end 



function vmecio_plot_vmec_comparison(
    equil, 
    vmec_data; 
    Nr = 10, 
    Nt = 8, 
    kwargs...
    )
    
    kwargs_dict = Dict(pairs(kwargs)) 

    py"""
    import desc 
    import desc.vmec

    def vmec():
        return desc.vmec.VMECIO.plot_vmec_comparison(
           $equil, 
           $vmec_data, 
           Nr = $Nr, 
           Nt = $Nt, 
           **$kwargs_dict
        )
    """
    output = py"vmec"()
end 


function vmecio_read_vmec_output(
    fname
) 
    py"""
    import desc 
    import desc.vmec

    def vmec():
        return desc.vmec.VMECIO.vmecio_read_vmec_output(
           $fname
        )
    """
    output = py"vmec"()
end 


function vmecio_save(
    eq, 
    path; 
    surfs = 128, 
    verbose = 1
) 
    py"""
    import desc 
    import desc.vmec

    def vmec():
        return desc.vmec.VMECIO.save(
           $eq, 
           $path, 
           surfs = $surfs, 
           verbose = $verbose
        )
    """
    output = py"vmec"()
end 

function vmecio_interpolate(
    Cmn, 
    Smn, 
    xm, 
    xn, 
    theta,   
    phi, 
    s = nothing, 
    si = nothing, 
    sym = true
    )
    
    py"""
    import desc 
    import desc.vmec

    if isinstance($Cmn, np.ndarray):
        new_Cmn = np.ascontiguousarray($Cmn)
        assert new_Cmn.flags['C_CONTIGUOUS']
    else:
        new_Cmn = $Cmn
    if isinstance($Smn, np.ndarray):
        new_Smn = np.ascontiguousarray($Smn)
        assert new_Smn.flags['C_CONTIGUOUS']
    else:
        new_Smn = $Smn

    if isinstance($xm, np.ndarray):
        new_xm = np.ascontiguousarray($xm)
        assert new_xm.flags['C_CONTIGUOUS']
    else:
        new_xm = $xm

    if isinstance($xn, np.ndarray):
        new_xn = np.ascontiguousarray($xn)
        assert new_xn.flags['C_CONTIGUOUS']
    else:
        new_xn = $xn

    if isinstance($theta, np.ndarray):
        new_theta = np.ascontiguousarray($theta)
        assert new_theta.flags['C_CONTIGUOUS']
    else:
        new_theta = $theta

    if isinstance($phi, np.ndarray):
        new_phi = np.ascontiguousarray($phi)
        assert new_phi.flags['C_CONTIGUOUS']
    else:
        new_phi = $phi

    if isinstance($si, np.ndarray):
        new_si = np.ascontiguousarray($si)
        assert new_si.flags['C_CONTIGUOUS']
    else:
        new_si = $si
        
    def vmec():
        return desc.vmec.VMECIO.vmec_interpolate(
           new_Cmn, 
           new_Smn, 
           new_xm, 
           new_xn, 
           new_theta, 
           new_phi, 
           s = new_s, 
           si = new_si, 
           sym = $sym 
        )
    """
    output = py"vmec"()
end 
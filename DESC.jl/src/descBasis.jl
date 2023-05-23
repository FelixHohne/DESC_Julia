using PyCall 

#=
Functions implementing desc.basis API
=#

function PowerSeries(
    L;
    sym = "even", 
) 
    py"""
    import desc 
    import desc.basis
    import numpy as np  

    def create_power_series_basis():
        return desc.basis.PowerSeries(
           $L, 
           sym = $sym
        )
    """
    output = py"create_power_series_basis"()
end 

function FourierSeries(
    N;
    NFP = 1, 
    sym = false
) 
    py"""
    import desc 
    import desc.basis
    import numpy as np  

    def create_fourier_series_basis():
        return desc.basis.FourierSeries(
            $N,
            NFP = $NFP, 
            sym = $sym
        )
    """
    output = py"create_fourier_series_basis"()
end 

function DoubleFourierSeries(
   M, 
   N; 
   NFP = 1, 
   sym = false
) 
    py"""
    import desc 
    import desc.basis
    import numpy as np  

    def create_double_fourier_series_basis():
        return desc.basis.DoubleFourierSeries(
            $M, 
            $N, 
            NFP = $NFP, 
            sym = $sym
        )
    """
    output = py"create_double_fourier_series_basis"()
end 

function ZernikePolynomial(
    L, 
    M;
    sym = false, 
    spectral_indexing = "ansi"
 ) 
     py"""
     import desc 
     import desc.basis
     import numpy as np  
 
     def create_zernike_polynomial_basis():
         return desc.basis.ZernikePolynomial(
             $L, 
             $M, 
             sym = $sym, 
             spectral_indexing = $spectral_indexing
         )
     """
     output = py"create_zernike_polynomial_basis"()
 end 
 
 function FourierZernikeBasis(
    L, 
    M,
    N; 
    NFP = 1, 
    sym = false, 
    spectral_indexing = "ansi"
 ) 
     py"""
     import desc 
     import desc.basis
     import numpy as np  
 
     def create_basis():
         return desc.basis.FourierZernikeBasis(
             $L, 
             $M,
             $N,  
             sym = $sym, 
             spectral_indexing = $spectral_indexing
         )
     """
     output = py"create_basis"()
 end 
 
function basis_evaluate(
    basis, 
    nodes; 
    derivatives = [0, 0, 0], 
    modes = nothing, 
    unique = false)
    py"""
    import desc 
    import desc.basis
    import numpy as np  
    def ev():
        if isinstance($nodes, np.ndarray):
            new_nodes = np.ascontiguousarray($nodes)
            assert new_nodes.flags['C_CONTIGUOUS']
        else:
            new_nodes = $nodes
        
        if isinstance($derivatives, np.ndarray):
            new_derivatives = np.ascontiguousarray($derivatives)
            assert new_derivatives.flags['C_CONTIGUOUS']
        else:
            new_derivatives = $derivatives
        
        if isinstance($modes, np.ndarray):
            new_modes = np.ascontiguousarray($modes)
            assert new_modes.flags['C_CONTIGUOUS']
        else:
            new_modes = $modes

        return $basis.evaluate(
            new_nodes, 
            derivatives=new_derivatives, 
            modes = new_modes, 
            unique = $unique
        ) 
    """
    output = py"ev"()
end

function basis_get_idx(
    basis;
    L = 0, 
    M = 0, 
    N = 0, 
    error = true
)
    py"""
    import desc 
    import desc.basis 
    def compute():
        $basis.get_idx(
            L = $L, 
            M = $L, 
            N = $N, 
            error = $error
        ) 
    """
    output = py"compute"()
end


# basis type one of `{PowerSeries, FourierSeries, DoubleFourierSeries, ZernikePolynomial, FourierZernikeBasis}`
function basis_load(
   basis_type :: String,
   load_from; 
   file_format = nothing
)
    py"""
    import desc 
    import desc.basis 
    def load():
        if $basis_type == "PowerSeries"
            return desc.basis.PowerSeries.load(
                $load_from, 
                file_format = $file_format
            ) 
        elif $basis_type == "FourierSeries":
            return desc.basis.FourierSeries.load(
                $load_from, 
                file_format = $file_format
            ) 
        elif $basis_type == "DoubleFourierSeries":
            return desc.basis.DoubleFourierSeries.load(
                $load_from, 
                file_format = $file_format
            ) 
        elif $basis_type == "ZernikePolynomial":
            return desc.basis.ZernikePolynomial.load(
                $load_from, 
                file_format = $file_format
            ) 
        elif $basis_type == "FourierZernikeBasis":
            return desc.basis.FourierZernikeBasis.load(
                $load_from, 
                file_format = $file_format
            ) 
    """
    output = py"load"()
end





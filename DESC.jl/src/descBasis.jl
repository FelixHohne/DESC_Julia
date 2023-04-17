using PyCall 

function jl_basis_power_series(
    L;
    sym = "even", 
) 
    py"""
    import desc 
    import desc.basis
    import numpy as np  

    def create_power_series_basis():
        return desc.basis.PowerSeries(
           L = $L, 
           even = $even
        )
    """
    output = py"create_power_series_basis"()
end 

function jl_basis_fourier_series(
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
            N = $N;
            NFP = $NFP, 
            sym = $sym
        )
    """
    output = py"create_fourier_series_basis"()
end 

function jl_basis_double_fourier_series(
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
            M = $M, 
            N = $N, 
            NFP = $NFP, 
            sym = $sym
        )
    """
    output = py"create_double_fourier_series_basis"()
end 

function jl_basis_zernike_polynomial(
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
             L = $L, 
             M = $M, 
             sym = $sym, 
             spectral_indexing = $spectral_indexing
         )
     """
     output = py"create_zernike_polynomial_basis"()
 end 
 
 function jl_basis_fourier_zernike_polynomial(
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
 
     def create_fourier_zernike_polynomial_basis():
         return desc.basis.FourierZernikePolynomial(
             L = $L, 
             M = $M,
             N = $N,  
             sym = $sym, 
             spectral_indexing = $spectral_indexing
         )
     """
     output = py"create_fourier_zernike_polynomial_basis"()
 end 
 
function change_basis(basis, N; NFP=nothing) 

    py"""
    import desc 
    import desc.basis 
    def change_resolution():
        basis.change_resolution(N = $N, NFP = $NFP) 
    """
    output = py"change_resolution"()
end 
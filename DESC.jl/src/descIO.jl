using PyCall

#=
Functios implement desc.io API
=#
function InputReader(;
    cl_args = nothing
)

    py"""
    import numpy as np
    import desc

    def input_reader():
        return desc.io.InputReader(
            cl_args = $cl_args
        )
    """
    output = py"input_reader"()
end


function io_parse_vmec_inputs(
    vmec_fname; 
    threshold = 0
)

py"""
import numpy as np
import desc

def parse_vmec_inputs():
    return desc.io.parse_vmec_inputs(
        $vmec_fname, 
        threshold = $threshold
    )
"""
output = py"parse_vmec_inputs"()
end

function io_vmec_to_desc_input(
    vmec_fname, 
    desc_fname
)

py"""
import numpy as np
import desc

def vmec_to_desc_input():
    return desc.io.vmec_to_desc_input(
        $vmec_fname, 
        $desc_fname
    )
"""
output = py"vmec_to_desc_input"()
end


function io_write_desc_input(
    filename,
    inputs;
    header = ""
)

py"""
import numpy as np
import desc

def write_desc_input():
    return desc.io.write_desc_input(
       $filename, 
       $inputs, 
       header = $header
    )
"""
output = py"write_desc_input"()
end
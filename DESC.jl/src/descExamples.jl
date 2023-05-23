using PyCall 

#=
Functions implement desc.examples API
=#
function examples_get(name; data = nothing) 
    py"""
    import desc
    import desc.examples
    def get():
        return desc.examples.get(
            $name, 
            data = $data
        )
    """
    output = py"get"()
end 


function examples_listall() 
    py"""
    import desc
    import desc.examples
    def listall():
        return desc.examples.listall()
    """
    output = py"listall"()
end 

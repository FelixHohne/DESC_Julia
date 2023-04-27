using PyCall 

function examples_get(name; data = nothing) 
    py"""
    import desc
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
    def listall():
        return desc.examples.listall()
    """
    output = py"listall"()
end 

using YAML

"""
Convert a nested Dict (from YAML) to the Hamster block format.
"""
function write_blocks(io::IO, data::Dict)
    for (section, content) in data
        println(io, "begin $section")
        if content isa Dict
            for (k, v) in content
                if v isa Vector
                    println(io, "  $k = ", join(v, " "))
                else
                    println(io, "  $k = $v")
                end
            end
        end
        println(io, "end\n")
    end
end

function yaml_to_blockfile(yaml_file::String, out_file::String)
    data = YAML.load_file(yaml_file)
    open(out_file, "w") do io
        write_blocks(io, data)
    end
end



path_yaml = ARGS[1]
path_hconf = ARGS[2]

yaml_to_blockfile(path_yaml, path_hconf)

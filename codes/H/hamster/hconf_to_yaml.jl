import YAML
import Hamster

function to_dict(c::Hamster.Config)
    # start with the options at top-level
    d = Dict("Options" => c.options)
    # merge each block into top-level
    for (k, v) in c.blocks
        d[k] = v
    end
    return d
end

path_hconf = ARGS[1]
path_yaml = ARGS[2]

config_hconf = Hamster.get_config(filename=path_hconf)

YAML.write_file(path_yaml, to_dict(config_hconf))
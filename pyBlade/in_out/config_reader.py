import configparser


def read(cfg_file):

    ## Get the config object ##
    cfg = configparser.ConfigParser(comment_prefixes="#", empty_lines_in_values=False)
    
    ### Read the parameters ###
    cfg.read(cfg_file)

    return cfg
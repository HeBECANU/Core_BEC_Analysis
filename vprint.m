function str = vprint(str,verbose,level)
    if verbose>level
        str = cli_header(str,level);
    end

end
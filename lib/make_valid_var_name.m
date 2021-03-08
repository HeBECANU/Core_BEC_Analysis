function name=make_valid_var_name(name)
% because
% matlab.lang.makeValidName('sadasd sd','ReplacementStyle','underscore') is broken
% i made my own
    regex_cmd='[^a-zA-Z0-9]{1}';
    name=regexprep(name,regex_cmd,'_');       
end

  
function str = num2eng(input, useSI, useMu, spaceAfterNum, sigFigs, fullName, useMinus, useInf, trailZeros, pad, cellOutput)
%============================================================================================================================================================
%   string = num2eng(input,[useSI],[useMu],[spaceAfterNum],[sigFigs],[fullName],...
%                          [useMinus],[useInf],[trailZeros],[pad],[cellOutput])
%   where input variables in square brackets are optional - i.e. The function accepts 
%   anywhere from 1 to 11 inputs, inclusive. Converts an input numerical value into 
%   an engineering-formatted string, in either scientific format with multiples-of-
%   three exponent, or using SI prefixes e.g. k, M, n, p etc.
%   OR changes the tick labels in a set of axes or a colorbar, or all axes in a
%   figure, to use engineering notation or SI units.
%
%   Alternative call syntax:
%   string = num2eng(number,optionsStruct), whereby the control options are passed in
%   a structure. Addional options available with this call syntax are: 
%   noExp4Tenths, expWidth, expSign, expForce, unit, noPlural, capE, smallCapE, 
%   fracMant, and axes. 
%   This syntax is explained in more detail below, starting on line 113.
%
%   Version: 5.2 (version history at end of file)
%
%   Please note that throughout this file, the term "string" is referring to a vector/array of characters, rather than MATLAB's "string" data type
%
%   Input variables: (NOTE: for logical inputs, it is acceptable to use 0/1 instead of true/false)
%
%       input:          The numeric value(s) to be converted to a string. If this is the only input, the string returned will be in scientific format with
%                       multiples-of-three exponent. e.g. num2eng(0.001) returns '1e-3' (not including quotes).
%
%                       "input" can also be a handle to some axes, colorbar, MATLAB figure, or an array of any of these types. The tick labels of the axes, 
%                       colorbar, or tick labels of all axes in the figure(s) if "input" is a figure handle or array of figure handles, will be formatted to
%                       use engineering notation. Listeners will be created to update the tick labels as necessary, if the ticks change later. By default,
%                       only the x-axis labels are processed, but the option 'axes' gives control over this - see line 179 of this comment box.                       
%
%                       If "input" is empty (i.e. []) or not numeric, the returned string will be '' (empty char), or ' ' if spaceAfterNum=true.
%
%                       If "input" is a vector or matrix, num2eng will return a cell matrix with the same "shape" as the input matrix, with one string per
%                       cell, or, if option cellOutput=false, a single string with each number separated by either a space or a new line, depending upon the
%                       shape of the input (see explanation of option cellOutput on line 99 for more info).
%
%                       Complex numbers are supported, but SI prefixes will not be used.
%
%       useSI:          Logical. When set to true, SI prefixes will be used for any number with magnitude in the range [1 y, 1000 Y), instead of scientific
%                       notation. Defaults to false if not specified.
%                       Examples:
%                       num2eng(0.001  ,true) returns '1 m'
%                       num2eng(2.56e29,true) returns '256e+27'
%                       num2eng(1.5    ,true) returns '1.5'
%                       num2eng(-1.2e5 ,1)    returns '-120 k'
%
%       useMu:          Logical. If set to true, also forces useSI to true. Then, if magnitude of the input number lies in the range [1e-6, 1e-3), the
%                       lower-case Greek mu (Unicode U+03BC) will be used. Otherwise, 'u' will be used here instead. Defaults to false if not specified.
%                       Examples:
%                       num2eng(10e-6,true)      returns '10 u'
%                       num2eng(10e-6,true,true) returns '10 <micro>', where <micro> is the lower-case Greek mu (Unicode U+03BC)
%
%       spaceAfterNum:  Logical. When set to true, a space will be inserted after the numeric part of the output string, after the 'e' if the output string
%                       is in scientific format, and in all cases where useSI is true, even for numbers with magnitude in the range [1, 1000). Defaults to
%                       false. Examples ('x' means 'don't care'/'not relevant' in the below examples):
%                       num2eng(0.001  ,true ,x    ,x)    returns '1 m'
%                       num2eng(0.001  ,false,false,true) returns '1e-3 '
%                       num2eng(2.56e29,true ,x    ,true) returns '256e+27 '
%                       num2eng(1.5    ,true ,x    ,true) returns '1.5 '
%
%       sigFigs:        Positive integer. Specifies how many significant figures should be in the output string. Defaults to 5 if not specified.
%
%       fullName:       Logical. If fullName is true, the full name of the SI prefix will be used, rather than just a letter. If set to true, forces useSI to
%                       true. Overrides the "useMu" option if that has been set. Defaults to false if not specified. Examples:
%                       num2eng(0.001       ,true,false,false,5,true) returns '1 milli'
%                       num2eng(2.56e29     ,true,false,false,5,true) returns '256e+27'
%                       num2eng(1.5         ,true,false,false,5,true) returns '1.5'
%                       num2eng(4.5645667e-6,1   ,1    ,1    ,6,1)    returns '4.56457 micro'
%
%       useMinus:       Logical. If useMinus is true, and the input number and/or the exponent is negative, the returned string will use the "proper" minus
%                       character (U+2212) instead of the default ascii hyphen-minus (U+002D). Defaults to false.
%                       Examples (in the below, <minus> is the Unicode character U+2212):
%                       num2eng(-0.001     ,true ,false,false,5,true)       returns '-1 milli'
%                       num2eng(-0.001     ,true ,false,false,5,true ,true) returns '<minus>1 milli'
%                       num2eng(-0.001     ,false,false,false,5,false,true) returns '<minus>1e<minus>3'
%                       num2eng(-0.001     ,false,false,false,5,false)      returns '-1e-3'
%                       num2eng(1000-0.001j,false,false,false,5,false)      returns '1e+3 - 1e-3i'
%                       num2eng(1000-0.001j,false,false,false,5,false,true) returns '1e+3 <minus> 1e<minus>3i'
%
%       useInf:         Logical. Set to true to use the infinity symbol (U+221E) when the input number is infinite. Defaults to false if not specified.
%                       Examples:
%                       num2eng(1/0, x, x, false, 5, x, x, true) returns '+<inf>', where <inf> is the infinity symbol (U+221E)
%
%       trailZeros:     Logical. Set to true to include any trailing zeros required to pad the output to the specified number of significant figures.
%                       Defaults to false. Examples:
%                       num2eng(1       ,1,0,1,5,0,0,0) returns '1 ', whilst
%                       num2eng(1       ,1,0,1,5,0,0,1) returns '1.0000 '
%                       num2eng(-1.34e-7,1,0,1,5,0,0,1) returns '-134.00 n'
%
%       pad:            Integer. Pads the output string to the specified width, by addition of spaces. A positive number gives a right-justified string
%                       (padding spaces are inserted before the number) whilst a negative number gives a left-justified string (padding spaces are inserted
%                       after the number and any prefix). A value of zero does not insert any padding. Defaults to zero. Examples:
%                       num2eng(-1.34e-7,1,0,1,5,0,0,0,10)  returns '    -134 n'
%                       num2eng(-1.34e-7,1,0,1,5,0,0,0,-10) returns '-134 n    '
%                       num2eng(-1.34e-7,1,0,1,5,0,0,1,10)  returns ' -134.00 n'
%
%       cellOutput:     Logical. Only used if "input" is a vector or 2D matrix. If set to true, the output is a cell array with one cell per element
%                       of the input number. If false, the output is a single string, separated by spaces (normal spaces, not non-breaking spaces) for each
%                       new column in the input, and newlines for each new row in the input. Defaults to true if not specified.
%                       Examples:
%                       num2eng([-0.001 1e6 12e5],true, true, false, 4, false, false, false, 0, true) will return the cell array:
%                       {'-1 m', '1 M', '1.2 M'}
%                       whereas
%                       num2eng([-0.001 1e6 12e5],true, true, false, 4, false, false, false, 0, false) returns:
%                       '-1 m 1 M 1.2 M'
%                       and
%                       num2eng([-0.001 1e6 12e5; 3 56e32 6.78e7],false, true, false, 4, false, false, true, 0, false) returns:
%                       '-1.000 m 1.000 M 1.200 M
%                        3.000 5.600e+33 67.80 M'
%
%   Using an options structure instead of individual option inputs (see from line 131 for the additional options available):
%       When num2eng was first developed, the individual option-input approach was selected in order to make function hints as helpful as possible when
%       writing code. However, as the number of options has grown, the function call has become unwieldy, especially if you only want to set one of the later
%       options and leave the others at the default value, or when reading the code later. This is where using the alternative syntax comes in. You can pass
%       num2eng an options structure as the second input. This structure should have anywhere from one to 20 fields, named as per the options listed above,
%       or the ten additional options listed below (starting on line 131). You can use the function hints to remind you of the option names for naming of
%       your option struct fields, for the first ten options.
%       For example:
%
%       num2eng(18.9e-10, struct('fullName',true) ) returns:
%       '1.89 nano'
%       and
%       num2eng([12.345563e-10 89e5 -0.00003; 0.1 14.5e2 77e8], struct('sigFigs',7, 'cellOutput',false, 'pad',10) ) returns:
%       '1.234556e-9     8.9e+6     -30e-6
%            100e-3    1.45e+3     7.7e+9'
%
%       Additional options available with this call syntax:
%
%       noExp4Tenths:   Logical. When set to true, if the magnitude of the input number lies in the range [0.1, 1), no exponent or SI prefix is used.
%                       Defaults to false; over-ridden by expForce. examples:
%                       num2eng(0.2)   returns '200e-3'
%                       num2eng(0.2,1) returns '200 m'
%                       num2eng(0.2,struct('noExp4Tenths',1))           returns '0.2'
%                       num2eng(0.2,struct('useSI',1,'noExp4Tenths',1)) returns '0.2'
%
%       expWidth:       Integer. Sets the desired character width of the exponent (excluding sign character). The exponent will be padded with leading zeros
%                       in order to give the requested width, but it will never be truncated. Examples:
%                       num2eng(2e6,struct('expWidth',2))   returns '2e+06'
%                       num2eng(2e121,struct('expWidth',2)) returns '20e+120'
%
%       expSign:        Logical. Set to false to give a sign character in the exponent only for negative exponents. Examples:
%                       num2eng(1e4) returns '10e+3'
%                       num2eng(1e4,struct('expSign',false)) returns '10e3'
%
%       expForce:       Logical. Forces the exponent to always be present, even if it is zero. Overrides the noExp4Tenths option. Examples:
%                       num2eng(1) returns '1'
%                       num2eng(1,struct('expForce',true)) returns '1e+0'
%                       num2eng(0.2,struct('noExp4Tenths',true,'expForce',true)) returns '200e-3'
%                       num2eng(0.2,struct('useSI',true,'noExp4Tenths',true,'expForce',true)) returns '200 m'
%
%       unit:           Character vector or cell array of character vectors. Provides a string to append to the end of the output string, typically a unit of
%                       measure such as "volt", "meter", "A", etc. When processing vector/array/matrix inputs, the provided unit will be appended to all
%                       output strings. Words are automatically pluralised if the output number is not unity (but this can be overridden if desired using the
%                       "noPlural" option - see below). For vector/array/matrix number inputs, you can also provide a cell array of strings with the same
%                       "shape" as the number input, if different units for different inputs are required. Examples:
%
%                       num2eng(0.08, struct('useSI',true,'unit','A'))      returns '80 mA'
%                       num2eng(0.08, struct('fullName',true,'unit','amp')) returns '80 milliamps'
%                       num2eng([0.08,5e4], struct('fullName',true,'unit',{{'amp','watt'}})) returns the 1 x 2 cell array {'80 milliamps', '50 kilowatts'}
%
%                       NOTE: when using "struct" to create the input structure, and you wish to specify 'unit' as a cell array, you must use a double set of
%                       curly braces as shown above - this is standard MATLAB syntax.
%
%       noPlural:       Logical. Set to true to override the automatic pluralising of unit strings:
%                       num2eng(0.08,struct('fullName',true,'unit','amp','noPlural',true)) returns '80 milliamp'
%
%       capE:           Logical. Set to true to use a capital "E" instead of lower-case "e" for the exponent part (when not using SI prefixes).
%
%       smallCapE:      Logical. Set to true to use a small capital "E" (unicode U+1D07) instead of lower-case "e" for the exponent part (when not using SI 
%                       prefixes). If true, overrides the "capE" setting.
%
%       fracMant:       Logical. Set to true to force the mantissa (first numeric part) of the output string to lie in the range [0.001, 1), instead of the
%                       more usual [1, 1000). This would typically be used when you want to be unambiguous regarding the accuracy of the stated output string
%                       i.e. 0.1 mm meaning 0.1 mm +/- 0.05 mm, whilst 0.10 mm means 0.1 mm +/- 0.005 mm. This would not be possible to convey under the
%                       standard notation, that would return 100 um. (see https://en.wikipedia.org/wiki/Engineering_notation#Overview for more).
%
%       axes:           Character vector. If "input" is a handle or array of handles to axes, or figures, this option sets which axes (x, y, or z) of the
%                       input axes are processed to format their tick labels. Simply list the axes you wish to process. e.g. set to 'y' to process only the
%                       y-axis, or 'yz' to process y-axis and z-axis. Also, 'off' is a valid setting - this reverts any previously-processed labels to the
%                       default and stops listening for axis tick changes.
%                       Examples:
%                       num2eng(gca,struct('axes','y'))               formats y-axis tick labels of current axes to use engineering notation
%                       num2eng(gca,struct('useSI',true,'axes','xy')) formats x and y axis of current axes to use numbers with SI prefixes
%                       num2eng(gcf,struct('axes','off'))             turns off any num2eng formatting on all axes and colorbars in the figure 
%
%
%   Output variables:
%       str:    The output string in engineering format, or cell-vector/matrix of strings if "input" is a vector/matrix. For vector/matrix input, see also 
%               the cellOutput option.
%
%   Written by Harry Dymond, Electrical Energy Management Group, University of Bristol, UK; harry.dymond@bris.ac.uk. If you find any bugs, please email me!
%   Developed with MATLAB R2018a running on Windows 10 and macOS Mojave. Inspired by a post by Jan Simon in the thread:
%   https://uk.mathworks.com/matlabcentral/answers/892-engineering-notation-printed-into-files
%   with further inspiration from Stephen Cobeldick's num2sip: 
%   https://uk.mathworks.com/matlabcentral/fileexchange/33174-number-to-scientific-prefix
%   and Roman Müller-Hainbach's NUM2ENG:
%   https://uk.mathworks.com/matlabcentral/fileexchange/63928-num2eng-fast-number-to-engineering-notation-si-conversion
%
%   The code herein may be freely shared and modified as long as this comment box is included in its entirety and any code changes and respective authors are
%   summarised here.
%
%============================================================================================================================================================

    %% Constants
    NL             = char(10);   %#ok<CHARTEN>, for backwards compatibility
    NO_BREAK_SPACE = char(160);
    %% Check inputs
    nInputs = nargin;
    if nInputs > 1 && isstruct(useSI) && isfield(useSI,'Processed')
        opts = useSI;
    else
        assert(nInputs>=1,['ERROR: At least one input argument is required. Required function syntax is:\n'...
                           '       outputString = '...
                           'num2eng(input,[useSI],[useMu],[spaceAfterNum],[sigFigs],[fullName],[useMinus],'...
                                         '[useInf],[trailZeros],[pad],[cellOutput],[expWidth],[expSign],[forceExp]);\n'...
                           '       where input variables in square brackets are optional\n'...             
                           'OR:    outputString = num2eng(input, optionStructure);\n' ...
                          ],[]);
                    %     1       2          3             4          5         6         7          8         9        10            11           12
        optionNames = {'useSI','useMu','spaceAfterNum','sigFigs','fullName','useMinus','useInf','trailZeros','pad','cellOutput','noExp4Tenths','expWidth'};
        defaults    = { false , false ,     false     ,    5    ,   false  ,  false   , false  ,   false    ,  0  ,    true    ,     false    ,     0    };
                    %                  13         14       15       16       17        18          19      20
        optionNames = [optionNames {'expSign','expForce','unit','noPlural','capE','smallCapE','fracMant','axes'}];
        defaults    = [defaults    {  true   ,  false   ,  ''  ,   false  , false,    false  ,   false  ,  'x' }];
                    %      1         2         3        4         5         6         7         8        9        10         11       12        13  
        type        = {'logical','logical','logical','scalar','logical','logical','logical','logical','scalar','logical','logical','scalar','logical',...
                       'logical','char|cell','logical','logical','logical','logical','char'};
                    %      14        15          16        17        18        19      20
        
        for i = 1:numel(optionNames), opts.(optionNames{i}) = defaults{i}; end
        if nInputs > 1
            if isstruct(useSI)
                assert(nInputs==2,'\nERROR: When using a structure to define the options, call syntax is:\n       num2eng(input, optionsStruct);',[]);
                inputOptions      = useSI;
                inputOptionFields = fieldnames(inputOptions);
                checkRange        = [];
                for i = 1:length(inputOptionFields)
                    [valid, idx]  = ismember(inputOptionFields{i}, optionNames);
                    if valid
                        opts.(inputOptionFields{i}) = inputOptions.(inputOptionFields{i});
                        checkRange(end+1) = idx; %#ok<AGROW>
                    else
                        fprintf(1,'[\bWARNING: Ignoring unknown option "%s"]\b\n',inputOptionFields{i});
                        % Undocumented [\b ... ]\b switch to print in orange, see: https://undocumentedmatlab.com/blog/another-command-window-text-color-hack
                    end
                end
            else
                    varIdx = 0;  set_opt(useSI);
                if nInputs > 2 , set_opt(useMu);         end
                if nInputs > 3 , set_opt(spaceAfterNum); end
                if nInputs > 4 , set_opt(sigFigs);       end
                if nInputs > 5 , set_opt(fullName);      end
                if nInputs > 6 , set_opt(useMinus);      end
                if nInputs > 7 , set_opt(useInf);        end
                if nInputs > 8 , set_opt(trailZeros);    end
                if nInputs > 9 , set_opt(pad);           end
                if nInputs > 10, set_opt(cellOutput);    end
                checkRange = 1:nInputs-1;
            end
            for varIdx = checkRange, opts.(optionNames{varIdx}) = check(opts.(optionNames{varIdx}),type{varIdx}); end %#ok<FXUP>
        end
        %% Deal with axes, figure, or colorbar handle inputs 
        if isa(input,'matlab.graphics.axis.Axes') || isa(input,'matlab.ui.Figure') || isa(input,'matlab.graphics.illustration.ColorBar')
            opts.Processed  = true;
            opts.cellOutput = true;
            opts.axes       = upper(opts.axes);
            barH            = matlab.graphics.illustration.ColorBar.empty;
            axesH           = matlab.graphics.axis.Axes.empty;
            if     isa(input,'matlab.graphics.axis.Axes')            , figsH = false; axesH = input; %#ok<ALIGN>
            elseif isa(input,'matlab.graphics.illustration.ColorBar'), figsH = false; barH  = input;    
            else                                                     , figsH = input;                end
            for figH = figsH(:)'
                if ~islogical(figH)
                    axesH = figH.Children(arrayfun(@(x)isa(x,'matlab.graphics.axis.Axes'),figH.Children)); 
                    barH  = figH.Children(arrayfun(@(x)isa(x,'matlab.graphics.illustration.ColorBar'),figH.Children));
                end
                for objH = [axesH(:)' barH(:)']
                    if isa(objH,'matlab.graphics.axis.Axes'), AXES = 'XYZ';
                    else                                    , AXES = '_';   end
                    for axis = AXES
                        if axis == '_', barChar = 's'; axis = ''; %#ok<FXSET>
                        else          , barChar = ''; end
                        if strcmp(opts.axes, 'OFF')
                            listener = getappdata(objH,['num2eng' axis 'Listener']);
                            if ~isempty(listener)
                                delete(listener)
                                setappdata(objH,['num2eng' axis 'Listener'],[])
                                objH.([axis 'TickLabel' barChar 'Mode']) = 'auto';
                            end
                        elseif contains(opts.axes, axis)
                            objH.([axis 'TickLabels']) = num2eng(objH.([axis 'Tick' barChar]), opts);
                            setappdata(objH,['num2engOld' axis 'Ticks'],objH.([axis 'Tick' barChar]));
                            listener = getappdata(objH,['num2eng' axis 'Listener']);
                            if isempty(listener)
                                listener = addlistener(objH.([axis 'Ruler']), 'MarkedClean', @(ruler, eventdata)update_tick_labels(ruler, axis, opts));
                                setappdata(objH,['num2eng' axis 'Listener'], listener);
                            else
                                % make sure listener callback is up-to-date with latest options, in case user called num2eng again with new ones
                                listener.Callback = @(ruler, eventdata)update_tick_labels(ruler, axis, opts);
                            end
                        end
                    end
                end
            end
            return
        end
    end

    %% Initialise space character
    if ~isempty(opts.unit), opts.spaceAfterNum = true; end
    if opts.spaceAfterNum, spaceChar = NO_BREAK_SPACE;
    else                 , spaceChar = '';  end
    %% Deal with edge cases
    if isempty(input) || ~isnumeric(input) || iscell(input)
        str = pad_str(spaceChar, opts.pad);
        return;
    end
    %% Check sig figs and exp width settings are valid
    if opts.sigFigs <= 0
        fprintf(1,'[\bWARNING: option "sigFigs" should be greater than zero.\nReverting to default value of %i]\b\n', defaults{4});
        opts.sigFigs = defaults{4};
    end
    if opts.useMu || opts.fullName, opts.useSI = true; end
    if opts.expWidth < 0
        fprintf(1,'[\bWARNING: option "expWidth" should be greater than or equal to zero.\nReverting to default value of %i]\b\n', defaults{12});
        opts.expWidth = defaults{12};
    end
    %% Make sure option "unit" is the right "shape"
    shortfall = 0;
    if ~isempty(opts.unit) && iscell(opts.unit) && numel(opts.unit)~=numel(input)
        shortfall = numel(input)-numel(opts.unit);
        if shortfall < 1 , issueStr  = 'Too many';                                                 %#ok<ALIGN>
                           actionStr = 'Deleting the last ';
        else             , issueStr  = 'Not enough'; 
                           actionStr = 'Using the last valid entry in "unit" for the final ';      end
        if abs(shortfall) == 1, actionStr  = [actionStr 'unit string.'];
        else                  , actionStr  = [actionStr num2str(abs(shortfall)) ' unit strings.']; end
        if numel(input) == 1 , elementStr = 'element';
        else                  , elementStr = 'elements';                                           end
        fprintf(1,'[\bWARNING: %s strings in option "unit" (input number has %i %s but unit has %i).\n%s\n]\b\n', ...
                  issueStr, numel(input), elementStr, numel(opts.unit), actionStr);
    elseif ~iscell(opts.unit)
        opts.unit = {opts.unit};
        shortfall = numel(input)-numel(opts.unit);
    end
    if shortfall
        opts.unit = opts.unit(:);
        if shortfall>0, opts.unit(end+1:end+shortfall) = repmat(opts.unit(end),1,shortfall);
        else          , opts.unit(end+shortfall+1:end) = [];                                    end 
    end
    if ~isequal(size(opts.unit), size(input)), opts.unit = reshape( opts.unit, size(input) ); end
    %% Initialise minus character
    if opts.useMinus, minusChar = char(8722); % "true" minus character U+2212
    else            , minusChar = '-';        end    
    %% Initialise exponent character
    if     opts.smallCapE , expChar = char(7431); %#ok<ALIGN> small capital E
    elseif opts.capE      , expChar = 'E';
    else                  , expChar = 'e';        end    
    %% Initialise SI prefix, if requested
    if opts.useSI
        expValue = [  24,   21,   18,   15,   12,    9,    6,    3,  0,   -3,   -6,   -9,  -12,  -15,  -18,  -21,  -24];
        if opts.fullName
            expName = {'yotta', 'zetta', 'exa', 'peta', 'tera', 'giga', 'mega', 'kilo', '', ...
                       'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', 'yocto'};
        else
            expName = {'Y', 'Z', 'E', 'P', 'T', 'G', 'M', 'k', '', 'm', 'u', 'n', 'p', 'f', 'a', 'z', 'y'};
            if opts.useMu, expName{11} = char(956); end % char(956) is the lower case Greek "mu" character
        end
        for i=1:17, expName{i} = [NO_BREAK_SPACE expName{i}]; end
        if ~opts.spaceAfterNum, expName{9} = ''; end
    else
        expValue = [];
        expName  = [];
    end
    %% Initialise trailing zeros sprintf format specifier
    if opts.trailZeros, trz = '#';
    else              , trz = '';  end
    %% Initialise infinity string
    if opts.useInf, infStr = [char(8734) spaceChar]; % infinity symbol
    else          , infStr = ['Inf'      spaceChar]; end
    %% Do the conversion
    if isscalar(input)
        str = num2eng_part1(input, opts.unit, opts, spaceChar, minusChar, expChar, expName, expValue, trz);
    elseif iscolumn(input) || isrow(input)
        % faster than using arrayfun
        if opts.cellOutput
            str{numel(input)} = '';
            for i=1:numel(input)
                str{i} = num2eng_part1(input(i), opts.unit(i), opts, spaceChar, minusChar, expChar, expName, expValue, trz);
            end
            if iscolumn(input), str = str'; end
        else
            str = num2eng_part1(input(1),opts.unit(1), opts, spaceChar, minusChar, expChar, expName, expValue, trz);
            if isrow(input), delim = ' ';
            else            , delim = NL; end
            for i = 2:numel(input)
                str = [str delim num2eng_part1(input(i),opts.unit(i), opts, spaceChar, minusChar, expChar, expName, expValue, trz)]; %#ok<AGROW>
            end
        end
    else
        if ndims(input)>2, opts.cellOutput = true; end  %#ok<ISMAT> <--- bug in mlint? ndims(number) can be == 2, and also ismat(number)==true
        if opts.cellOutput
            % fastest implementation for array input would be to loop explicitly, but arrayfun is much easier and matrix inputs are the least likely
            str = arrayfun(@(x,u)num2eng_part1(x, u, opts, spaceChar, minusChar, expChar, expName, expValue, trz),input,opts.unit,'UniformOutput',false);
        else
            nRows = size(input,1);
            nCols = size(input,2);
            str = ''; delim = '';
            for row = 1:nRows
                for col = 1:nCols
                    str   = [str delim ...
                             num2eng_part1(input(row,col),opts.unit(row,col), opts, spaceChar, minusChar, expChar, expName, expValue, trz)]; %#ok<AGROW>
                    delim = ' ';
                end
                if row<nRows, str = [str NL]; end %#ok<AGROW>
                delim = '';
            end
        end
    end

    %% String-forming subroutines
    % many subroutines do not make use of variables shared with the calling function, for performance reasons - see:
    % https://uk.mathworks.com/matlabcentral/answers/99537-which-type-of-function-call-provides-better-performance-in-matlab
    function str = num2eng_part1(input, unit, opts, spaceChar, minusChar, expChar, expName, expValue, trz)
        if ~isreal(input)
            %%     Handle imaginary numbers
            opts.useSI = false;
            str        = num2eng_part2(real(input),'', opts, '', minusChar, expChar, expName, expValue, trz);
            if imag(input) > 0, signChar = '+';
            else               , signChar = minusChar; input = -1*input; end
            str = [str NO_BREAK_SPACE signChar NO_BREAK_SPACE ...
                   num2eng_part2(imag(input),'', opts, spaceChar, minusChar, expChar, expName, expValue, trz) 'i' spaceChar unit{1}];
            return
        else
            %%     Handle all other cases
            str = num2eng_part2(input, unit{1}, opts, spaceChar, minusChar, expChar, expName, expValue, trz);
        end
    end

    function str = num2eng_part2(number, unit, opts, spaceChar, minusChar, expChar, expName, expValue, trz)
        %%     Handle special cases
        if     isempty(number) || ~isnumeric(number) || iscell(number), str = pad_str(spaceChar              , opts.pad); return
        elseif number==0                                              , str = pad_str(['0'    spaceChar unit], opts.pad); return
        elseif isnan(number)                                          , str = pad_str(['NaN'  spaceChar]     , opts.pad); return
        elseif isinf(number),                            if number>0  , str = pad_str(['+'       infStr unit], opts.pad); 
                                                         else         , str = pad_str([minusChar infStr unit], opts.pad); end, return
        end
        %%     Handle all other cases
        if opts.noExp4Tenths && abs(number)<1 && abs(number)>=0.1 && ~opts.expForce
            mantissa = num_to_n_sig_figs(number, opts.sigFigs);
            exponent = 0;
            expIndex = 0;
        else
            exponent = 3 * floor(log10(abs(number)) / 3);
            mantissa = num_to_n_sig_figs(number / (10 ^ exponent), opts.sigFigs);

            if abs(mantissa)>=1000
                exponent = exponent+3;
                mantissa = mantissa/1000;
            end
            
            if opts.fracMant
                exponent = exponent+3;
                mantissa = mantissa/1000;                
            end

            if opts.useSI, expIndex = (exponent == expValue); end
        end
        
        sigDigits   = max(opts.sigFigs, 1+floor(log10(abs(mantissa))));
        sprintfArgs = {opts.sigFigs, mantissa, spaceChar};
        if opts.trailZeros && opts.sigFigs <= 1+floor(log10(abs(mantissa))), sprintfFormat = '%.*g';
        else                                                               , sprintfFormat = ['%' trz '.*g']; end
        if opts.useSI && any(expIndex)
            sprintfFormat          = [sprintfFormat '%s'];
            sprintfArgs([1, 3, 4]) = {sigDigits, expName{expIndex}, ''};
        elseif exponent~=0 || opts.expForce
            sprintfArgs([1, 3, 4]) = {sigDigits, exponent, spaceChar};
            sprintfFormat          = [sprintfFormat expChar '%'];
            if opts.expSign, sprintfFormat(end+1) = '+'; end
            if opts.expWidth>0
                if opts.expSign || exponent<0, expWidth = opts.expWidth+1; 
                else                         , expWidth = opts.expWidth;  end
                sprintfFormat = [sprintfFormat '0#*'];
                sprintfArgs   = [sprintfArgs(1:2) expWidth sprintfArgs(3:4)];
            end
            sprintfFormat(end+1) = 'd';
        end
        
        str = (sprintf([sprintfFormat '%s' unit], sprintfArgs{:}));
        if length(unit)>1 && mantissa~=1 && unit(end)~='s' && ~opts.noPlural
            str(end+1) = 's';
        end
        if opts.useMinus
            if (number<0)  , str(1) = minusChar; end
            if (exponent<0), str    = strrep(str,'-',minusChar); end
        end
        str = pad_str(str, opts.pad);
        
        function num = num_to_n_sig_figs(num, n)
            num = round(num, n, 'significant');
        end
    end

    function str = pad_str(str, padLength)
        if      padLength > 0 && length(str)<   padLength, str = [repmat( ' ', 1,padLength-length(str) ) str];
        elseif  padLength < 0 && length(str)<-1*padLength, str = [str repmat( ' ', 1,-1*padLength-length(str) )]; end        
    end

    %% Input processing subroutines
    function set_opt(inputOption)
        varIdx = varIdx+1; opts.(optionNames{varIdx}) = inputOption;
    end

    function var = check(var, type)
        EXCEPTION_ERROR_START = '\nERROR: %s option (%s) should be a ';
        switch type
            case 'logical'
                if ~islogical(var)
                    try var = logical(var); catch, end
                end
                try if isscalar(var), return, end, catch, end
                errorMsgEnd = 'logical scalar\n';
            case 'scalar'
                if isnumeric(var) && isscalar(var) && isreal(var)
                    var = round(var); return;
                end
                errorMsgEnd = 'scalar integer\n';
            case {'char','char|cell'}
                if ischar(var), return, end
                errorMsgEnd = 'character vector';
                if contains(type,'cell')
                    if iscell(var)
                        areChars = cellfun(@(x)ischar(x),var);
                        if all(areChars(:)), return, end
                    end
                    errorMsgEnd = [errorMsgEnd ' or cell array of character vectors\n'];
                else
                    errorMsgEnd = [errorMsgEnd '\n'];
                end
        end
        throwAsCaller( MException('num2eng:inputError',[EXCEPTION_ERROR_START errorMsgEnd], num2ord(varIdx), optionNames{varIdx}) );
    end

    function str = num2ord(num)
        switch num
            case  1, str = 'First'      ;
            case  2, str = 'Second'     ;
            case  3, str = 'Third'      ;
            case  4, str = 'Fourth'     ;
            case  5, str = 'Fifth'      ;
            case  6, str = 'Sixth'      ;
            case  7, str = 'Seventh'    ;
            case  8, str = 'Eighth'     ;
            case  9, str = 'Ninth'      ;
            case 10, str = 'Tenth'      ;
            case 11, str = 'Eleventh'   ;
            case 12, str = 'Twelfth'    ;
            case 13, str = 'Thirteenth' ;
            case 14, str = 'Fourteenth' ;
            case 15, str = 'Fifteenth'  ;
            case 16, str = 'Sixteenth'  ;
            case 17, str = 'Seventeenth';
            case 18, str = 'Eighteenth' ;
            case 19, str = 'Nineteenth' ;
            case 20, str = 'Twentieth'  ;
        end
    end
end

%% Tick labels updater
function update_tick_labels(ruler, axis, opts)
    persistent busy
    if busy, return; end
    busy = true; %#ok<NASGU>
    objH = ruler.Parent;
    if isempty(axis), barChar = 's';
    else            , barChar = ''; end
    if isequal(objH.([axis 'Tick' barChar]), getappdata(objH,['num2engOld' axis 'Ticks'])), busy = false; return, end
    objH.([axis 'TickLabels']) = num2eng(objH.([axis 'Tick' barChar]), opts);
    setappdata(objH,['num2engOld' axis 'Ticks'],objH.([axis 'Tick' barChar]));
    busy = false;
end

%============================================================================================================================================================
%
%   Version history:
%
%   1.0         -- First release.
%   1.0.1       -- Documentation change.
%   1.0.2       -- Documentation change.
%   1.0.3       -- Option "fullName" now overrides option "useMu".
%               -- If input is a column vector, output will now be a column vector of cells instead of a row vector of cells.
%   1.0.4       -- No longer requires the Image Processing Toolbox.
%   1.0.6       -- No longer uses varargin, to make function hints more useful. Inputs 2 to 6 remain optional.
%   1.1         -- Added the option to use the true minus character (U+2212) in returned strings. The new option is added in such a way that the function
%                  will remain backwards compatible - code using earlier versions of num2eng does not need to be changed.
%   1.1.1       -- Icon update on FEX.
%   1.1.2       -- Minor code cleanup (no changes in functionality).
%   1.1.3       -- If input number is minus infinity, the useMinus option is now respected.
%   1.1.4       -- Reverted to using char(956) for micro; this is the Greek letter mu, whilst char(181) is the deprecated "micro" symbol, although in many
%                  fonts these code-points both use the same glyph.
%               -- Added commas to code as suggested by the new code parser in MATLAB 2018. This will result in MATLAB versions prior to this suggesting that
%                  the commas can/should be removed.
%   1.1.5       -- Extremely minor code cleanup/optimisation.
%   2.0         -- Added the "cellOutput" option, for use with vector/matrix inputs.
%               -- Added the "trailZeros" option.
%               -- Added the "pad" option.
%               -- Added option to use a structure to pass in the options, instead of passing options in as a list of separate variables.
%               -- Function remains compatible with any previous-working calls, so code already using function does not need to be changed.
%               -- If inputs useMu or fullName are set to true this now forces useSI to true.
%               -- Now uses non-breaking space between number and prefix.
%   2.1         -- Now uses repmat to generate space padding (no change in functionality; code just easier to read).
%   2.2         -- 06 May 2020
%               -- Infinite inputs now return strings where Inf is correctly capitalised.
%               -- Added the "useInf" option. This has been added as option 7, meaning options trailZeros, pad, cellOutput from versions 2.0 through 2.1 have
%                  "moved" - any code using these options in a syntax type 1 call (i.e. not using a settings struct) will have to be updated. Apologies for
%                  any inconvenience caused.
%   2.3         -- 06 May 2020
%               -- Changed mechanism for addition of optional space at end of string. Will have negligible impact on performance, just makes code a bit more
%                  compact.
%               -- Other minor code cleanups and documentation fixes.
%   3.0         -- 28 May 2020
%               -- Documentation fixes.
%               -- Fixed a bug when using the trailZeros option.
%               -- Option pad now respected in all cases.
%               -- Added options noExp4Tenths, expWidth, expSign, expForce, unit, and noPlural. (only available with optionStruct call syntax).
%   3.1         -- 30 May 2020
%               -- Minor code cleanup (no changes to functionality).
%   4.0         -- 04 June 2020
%               -- Added options capE, smallCapE, and fracMant.
%               -- Small performance (speed) enhancement.
%   4.1         -- 12 June 2020
%               -- Now automatically ensures that there is always a space between number and its units, if units have been specified.
%               -- Simplified input checking method (minor code cleanup).
%               -- 14 June 2020: Retrospectively changed version numbers (3.2->4.0, 3.3->4.1) as new features were added in what was originally called 3.2,
%                  so that version should have been called 4.0.
%   5           -- 03 July 2020 
%               -- Added ability to process axes tick labels, including listeners to update the tick labels if the ticks change after first call to num2eng
%   5.1         -- 03 July 2020
%               -- Added support for colorbars
%   5.2         -- 04 July 2020
%               -- Can now handle arrays of figure handles
%
%============================================================================================================================================================

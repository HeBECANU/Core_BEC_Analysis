function [f,df] = f2wl(varargin)
% Converts wavlength to frequency, and vice versa, with (absolute) uncertainty 
% The conversion is reciprocal, so this function converts both ways.
%   Inputs:     Wavelength (or frequency) in meters (or Hz)
%         (opt) Uncertainty (in same units)
%   Outputs:    Frequency (or wavelength) in Hz (or m)
%               Uncertainty (in same units) (iff also passed as input)
% TODO
% example usage
% test script

    c = 299792458; %Update as necessary
    f = c./varargin{1};
    if nargin>1
        df = c*varargin{2}./varargin{1}.^2;
    end
end
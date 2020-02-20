

%most of the test functonality comes from https://au.mathworks.com/matlabcentral/fileexchange/52832-num2sepstr
function tests = test_num_with_sep
    %test_num_with_sep
    tests = functiontests(localfunctions);
    timeit(@() num_with_sep(123456.012345,'format','%f'))
end


function testGeneric(testCase)
    verifyEqual(testCase,num_with_sep(1000),'1,000')
    verifyEqual(testCase,num_with_sep(1234.5),'1,234.5') % trims trailing zeros with auto format
    verifyEqual(testCase,num_with_sep(1234.5,'format','%f'),'1,234.500,000') % but not when format is specified
    verifyEqual(testCase,num_with_sep(1234.5,'format','%.2f'),'1,234.50')
    verifyEqual(testCase,num_with_sep(1234.5,'format','%.4f'),'1,234.500,0')
    verifyEqual(testCase,num_with_sep(123456789.5,'format','%.0f'),'123,456,790')
end

function testZero(testCase)
    verifyEqual(testCase,num_with_sep(0),'0')
end
function testInteger(testCase)
    verifyEqual(testCase,num_with_sep(1),'1')
    verifyEqual(testCase,num_with_sep(-1),'-1')
    verifyEqual(testCase,num_with_sep(100),'100')
    verifyEqual(testCase,num_with_sep(-100),'-100')
end
function test_trim_zeros_default_format(testCase)
    verifyEqual(testCase,num_with_sep(1),'1')
    verifyEqual(testCase,num_with_sep(1234),'1,234')
    verifyEqual(testCase,num_with_sep(1234.00001),'1,234')
    verifyEqual(testCase,num_with_sep(1234.0100),'1,234.01')
    verifyEqual(testCase,num_with_sep(1234.01001),'1,234.01')
end

function test_trim_zeros_custom_format(testCase)
    verifyEqual(testCase,num_with_sep(1.000123000456,...
                        'format','%.15f',...
                        'remove_trailing_zero',0,...
                        'add_sep_after_decimal',1),...
                '1.000,123,000,456,000')
    verifyEqual(testCase,num_with_sep(1.000123000456,...
                        'format','%.15f',...
                        'remove_trailing_zero',1,...
                        'add_sep_after_decimal',1),...
                '1.000,123,000,456')
end


function testFormatSpecfiers(testCase)
    verifyEqual(testCase,num_with_sep(1234.6,'format','%.4f'),'1,234.600,0')
    verifyEqual(testCase,num_with_sep(1234.6,'format','%.0f'),'1,235')
    verifyEqual(testCase,num_with_sep(1234.6,'format','%+.4f'),'+1,234.600,0')
end
function testLongNumbers(testCase)
    verifyEqual(testCase,num_with_sep(1234567890123456),'1,234,567,890,123,456')
    verifyEqual(testCase,num_with_sep(1234567890123456,'format','%.0f'),'1,234,567,890,123,456')
    verifyEqual(testCase,num_with_sep(12345678901234567890),'12,345,678,901,234,567,168')
    verifyEqual(testCase,num_with_sep(12345678901234,'format','%.0f'),'12,345,678,901,234')
end
function testScientificNotation(testCase)
    verifyEqual(testCase,num_with_sep(126,'format','%.0g'),num2str(126,'%.0g'))
    verifyEqual(testCase,num_with_sep(126,'format','%.5g'),num2str(126,'%.5g'))
    verifyEqual(testCase,num_with_sep(126,'format','%.0e'),num2str(126,'%.0e'))
    verifyEqual(testCase,num_with_sep(126,'format','%.5e'),'1.260,00e+02')
    verifyEqual(testCase,num_with_sep(1235674786,'format','%.0g'),num2str(1235674786,'%.0g'))
    verifyEqual(testCase,num_with_sep(123456.012345,'format','+%.14e'),'+1.234,560,123,450,00e+05')
end

function test_comma_after_zero(testCase)
verifyEqual(testCase,num_with_sep(123456.012345,'format','%f','add_sep_after_decimal',1),'123,456.012,345')
end






%formated_str=num_with_sep(123456789.123456789,1,'%.15f')






function tests = test_num_with_comma
    %test_num_with_comma
    tests = functiontests(localfunctions);
    timeit(@() num_with_comma(123456.012345,'%f',1))
end


function testGeneric(testCase)
    verifyEqual(testCase,num_with_comma(1000),'1,000')
    verifyEqual(testCase,num_with_comma(1234.5),'1,234.5') % trims trailing zeros with auto format
    verifyEqual(testCase,num_with_comma(1234.5,'%f'),'1,234.500000') % but not when format is specified
    verifyEqual(testCase,num_with_comma(1234.5,'%.2f'),'1,234.50')
    verifyEqual(testCase,num_with_comma(1234.5,'%.4f'),'1,234.5000')
    verifyEqual(testCase,num_with_comma(123456789.5,'%.0f'),'123,456,790')
end

function testZero(testCase)
    verifyEqual(testCase,num_with_comma(0),'0')
end
function testInteger(testCase)
    verifyEqual(testCase,num_with_comma(1),'1')
    verifyEqual(testCase,num_with_comma(-1),'-1')
    verifyEqual(testCase,num_with_comma(100),'100')
    verifyEqual(testCase,num_with_comma(-100),'-100')
end
function test_trim_zeros_default_format(testCase)
    verifyEqual(testCase,num_with_comma(1),'1')
    verifyEqual(testCase,num_with_comma(1234),'1,234')
    verifyEqual(testCase,num_with_comma(1234.00001),'1,234')
    verifyEqual(testCase,num_with_comma(1234.0100),'1,234.01')
    verifyEqual(testCase,num_with_comma(1234.01001),'1,234.01')
end

function test_trim_zeros_custom_format(testCase)
    verifyEqual(testCase,num_with_comma(1.000123000456,'%.15f',1,0),'1.000,123,000,456,000')
    verifyEqual(testCase,num_with_comma(1.000123000456,'%.15f',1,1),'1.000,123,000,456')
end


function testFormatSpecfiers(testCase)
    verifyEqual(testCase,num_with_comma(1234.6,'%.4f'),'1,234.6000')
    verifyEqual(testCase,num_with_comma(1234.6,'%.0f'),'1,235')
    verifyEqual(testCase,num_with_comma(1234.6,'%+.4f'),'+1,234.6000')
end
function testLongNumbers(testCase)
    verifyEqual(testCase,num_with_comma(1234567890123456),'1,234,567,890,123,456')
    verifyEqual(testCase,num_with_comma(1234567890123456,'%.0f'),'1,234,567,890,123,456')
    verifyEqual(testCase,num_with_comma(12345678901234567890),'12,345,678,901,234,567,168')
    verifyEqual(testCase,num_with_comma(12345678901234,'%.0f'),'12,345,678,901,234')
end
function testScientificNotation(testCase)
    verifyEqual(testCase,num_with_comma(126,'%.0g'),num2str(126,'%.0g'))
    verifyEqual(testCase,num_with_comma(126,'%.5g'),num2str(126,'%.5g'))
    verifyEqual(testCase,num_with_comma(126,'%.0e'),num2str(126,'%.0e'))
    verifyEqual(testCase,num_with_comma(126,'%.5e'),num2str(126,'%.5e'))
    verifyEqual(testCase,num_with_comma(1235674786,'%.0g'),num2str(1235674786,'%.0g'))
    verifyEqual(testCase,num_with_comma(123456.012345,'+%.14e',1),'+1.234,560,123,450,00e+05')
end

function test_comma_after_zero(testCase)
verifyEqual(testCase,num_with_comma(123456.012345,'%f',1),'123,456.012,345')
end






%formated_str=num_with_comma(123456789.123456789,1,'%.15f')


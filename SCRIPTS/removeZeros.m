function output = removeZeros(input, epsilon)
% Remove all entries from input with abs(entry) < epsilon and return 
% resulting vector
output = input(abs(input) > epsilon);

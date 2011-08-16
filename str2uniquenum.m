function unique = str2uniquenum(str)

% First, convert the string into an array of char values
str_a = [];
for i=1:length(str)
  str_a(i) = double(str(i));
end

% Pad it so this array has a length of 4n, n integer
str_a = [ str_a zeros(4*(ceil(length(str_a)/4) - length(str_a)/4),1) ];

% For convenience, reshape it into a 4*n matrix
str_a = reshape(str_a, length(str_a)/4, 4);

% Iterate through the matrix, converting each row into a unique number
str_b = [];
for i=1:size(str_a, 1)
  str_b(i) = str_a(i,1) * (128^3) ...
           + str_a(i,2) * (128^2) ...
           + str_a(i,3) * (128^1) ...
           + str_a(i,4);
end

unique = str_b;
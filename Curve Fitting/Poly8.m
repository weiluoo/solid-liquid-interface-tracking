function vals = Poly8(parameters, x)

vals = ones(size(x)) .* parameters(9);
for i = 8:-1:1
    vals = vals + parameters(9 - i) .* (x .^ i);
end

end
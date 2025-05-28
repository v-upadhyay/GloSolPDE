terminal_state = ones(10)
terminal_array = zeros(10,10)
error_array = zeros(10, 10)
for i in 1:10
    error_array[:, i] = terminal_array[:, i] - terminal_state
end
print(error_array)
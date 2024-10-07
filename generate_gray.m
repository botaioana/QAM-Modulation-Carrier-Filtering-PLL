function gray_code = generate_gray(n)
    if n == 1
        gray_code = {'0', '1'};
    else
        previous_gray = generate_gray(n-1);
        gray_code = [strcat('0', previous_gray), strcat('1', flip(previous_gray))];
    end
end
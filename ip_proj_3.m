% Read the input image
input_image = imread('Knee.pgm');

%  Declare the variables
D1 = 20; % First cutoff frequency
D2 = 160;  % Second cutoff frequency
[M, N] = size(input_image);
P = 2*M - 1;
Q = 2*N - 1;


% Zero-padding the input image and centering the low frequencies
padded_image = zeros(P, Q);
padded_image(1:M, 1:N) = double(input_image);

for x = 1:P
    for y = 1:Q
        padded_image(x, y) = padded_image(x, y) * (-1)^(x+y);
    end
end

% Compute the 2D DFT
F = dft2d(padded_image);

% Define the ideal low-pass filter for first cut off frequncy H1(u, v)
H1 = zeros(P, Q);
for u = 1:P
    for v = 1:Q
        D0 = sqrt((u - P/2)^2 + (v - Q/2)^2);
        if D0 <= D1
            H1(u, v) = 1;
        end
    end
end

% Define the ideal low-pass filter for second cut off frequency for H2(u, v)
H2 = zeros(P, Q);
for u = 1:P
    for v = 1:Q
        D0 = sqrt((u - P/2)^2 + (v - Q/2)^2);
        if D0 <= D2
            H2(u, v) = 1;
        end
    end
end

% Apply the low-pass filter for first cut off frequncy H1(u, v)
G1 = H1 .* F;

% Apply the low-pass filter for second cut off frequncy H2(u, v)
G2 = H2 .* F;

% Compute the inverse DFT for the first cutoff frequency H1(u, v)
i1 = idft2d(G1);

% Compute the inverse DFT for the second cutoff frequency H2(u, v)
i2 = idft2d(G2);

% Take the real part and crop to the original size for the first cutoff frequency H1(u, v)
image_output_1 = real(i1(1:M, 1:N));

% Take the real part and crop to the original size for the second cutoff frequency H2(u, v)
image_output_2 = real(i2(1:M, 1:N));

% Save the output images for both first and second cuttoff frequencies
imwrite(uint8(image_output_1), 'image_output_1.pgm');
imwrite(uint8(image_output_2), 'image_output_2.pgm');

% 2D DFT
function F = dft2d(z)
    [M, N] = size(z);
    F = zeros(M, N);
    for u = 1:M
        F(u, :) = dft(z(u, :));
    end
    for v = 1:N
        F(:, v) = dft(F(:, v).');
    end
end

% 2D inverse DFT
function z = idft2d(F)
    [M, N] = size(F);
    z = zeros(M, N);
    for m = 1:M
        z(m, :) = idft(F(m, :));
    end
    for n = 1:N
        z(:, n) = idft(z(:, n).');
    end
end

% 1D DFT
function Z = dft(x)
    N = length(x);
    Z = zeros(1, N);
    for k = 0:N-1
        for n = 0:N-1
            Z(k+1) = Z(k+1) + x(n+1) * exp(-1i * 2 * pi * k * n / N);
        end
    end
end

% 1D inverse DFT
function z = idft(X)
    N = length(X);
    z = zeros(1, N);
    for n = 0:N-1
        for k = 0:N-1
            z(n+1) = z(n+1) + X(k+1) * exp(1i * 2 * pi * k * n / N);
        end
        z(n+1) = z(n+1) / N;
    end
end

% pos_sphere = @(r, theta, phi) r.*[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

r = 1;
% theta = (0:180)*pi/180;
phi = (-180:180)*pi/180;

pos = @(r, phi) r.*[cos(phi); sin(phi); zeros(size(phi))];

pos_1 = pos(r,phi);


phi_source = 150*pi/180;
pos_source = pos(r, phi_source);
amplitude_source = 20;
figure
plot3(pos_1(1,:), pos_1(2,:), phi)
view(0,90)
hold on
plot3(pos_source(1), pos_source(2), phi_source, 'x')
grid on

dx = 0.15;
dy = 0.15;
Nx = 2;
Ny = 1;
A1 = array (0,0,0, 0, dx, dy, Nx, Ny);

plot_array(A1.A)

fs = 64000;
var_data = 3; % 1 - sine, 2 - noise, 3 - beepbeep
T = 10;
f = 1000;
y_raw = audio_data(var_data, fs, T, f);
y_A1 = A1.record(fs,y_raw,pos_source,amplitude_source);
[beta, Detected_powers] = A1.estimate_AOA(y_A1,fs);
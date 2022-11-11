% white Gaussian signal
y = randn(1024,1);
% OR real data
y = audioread('beepbeep.wav');
% sampling frequency in hertz
f=48000;
% microphone distance in meter
d = 0.175;
% speed of sound in meter/second
v = 343;
% maximum TDOA in samples
dn_max = round(f*d/v);
% array positions in meter
arrays = [2.69, 0.19; 
           5.19, 3.21;
           2.70, 5.27;
           0.37, 3.54];
% simulated recording
x=nktp_sim(f,y,[2, 3],ones(1,8),5*ones(1,8),0, arrays, d);
% size of the data
N = size(x,1);

[p_hat, theta] = estimate_position(x, dn_max, N, arrays);
%%
figure; hold on
scatter(p_hat(1), p_hat(2), 50, 'ok', 'LineWidth', 1);
scatter(2, 3,100, 'xr', 'LineWidth', 2);

cm = colormap(hsv(4));
l = 4;
for i = 1:size(arrays, 1)
    x0 = arrays(i,1);
    y0 = arrays(i,2);
    scatter(x0,y0,200, 'o', 'MarkerEdgeColor', cm(i,:), 'LineWidth', 2);

    x1 = x0 + l*cos(theta(i));
    y1 = y0 + l*sin(theta(i));
    plot([x0,x1],[y0,y1], 'color', cm(i,:), 'LineWidth', 2)
    
end
% xlabel('$x$','Interpreter','LaTex')
% ylabel('$y$','Interpreter','LaTex')
% matlab2tikz('position_sim.tikz');
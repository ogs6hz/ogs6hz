path(path,'/Users/Livvy/Documents/BME 8730/Field_II_ver_3_24_windows');
%path(path,'/home/user/field_II/m_files');
field_init(-1)
% Set initial parameters
% Define the transducer
rect=[1 0/1000 0/1000 0 2/1000 0/1000 0 2/1000 3/1000 0 0/1000 3/1000 0 ...
1 2/1000 3/1000 1/1000 1.5/1000 0
1 -2/1000 0/1000 0 0/1000 0/1000 0 0/1000 3/1000 0 -2/1000 3/1000 0 ...
1 2/1000 3/1000 -1/1000 1.5/1000 0
1 -2/1000 -3/1000 0 0/1000 -3/1000 0 0/1000 0/1000 0 -2/1000 0/1000 0 ...
1 2/1000 3/1000 -1/1000 -1.5/1000 0
1 0/1000 -3/1000 0 2/1000 -3/1000 0 2/1000 0/1000 0 0/1000 0/1000 0 ...
1 2/1000 3/1000 1/1000 -1.5/1000 0];
center=[0 0 0];
focus=[0 0 0]/1000;
Th = xdc_rectangles (rect, center, focus);
xdc_show(Th)
points = [2.8,3.75,5;2,2.25,5;1.2,2.25,5;1.2,0,5;0,0,5];
x_a = [1.4,1,0.6,0.6,0];
y_b = [1.25,.75,.75,0,0];
[h, start_time] = calc_h(Th,points);
length_h = size(h,1);
t = start_time;
for i = 1:(length_h-1)
    t = [t; start_time+ele_size*i];
end
figure;
for k = 1:5
    subplot(5,1,k)
    plot(t,h(:,k));
    title(['Impulse Response at z/a = 10, b/a = 1.5, x/a = ', num2str(x_a(k)), ' and y/b = ',num2str(y_b(k))])
    xlabel('t')
    ylabel('h(x,t)')
    ylim([0 1])
    hold on
end
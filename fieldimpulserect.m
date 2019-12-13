path(path,'/Users/Livvy/Documents/BME 8730/Field_II_ver_3_24_windows');
%path(path,'/home/user/field_II/m_files');
field_init(-1)
% Set initial parameters
% Define the transducer
% % rect=[1 0/1000 0/1000 0 2/1000 0/1000 0 2/1000 3/1000 0 0/1000 3/1000 0 ...
% % 1 2/1000 3/1000 1/1000 1.5/1000 0
% % 1 -2/1000 0/1000 0 0/1000 0/1000 0 0/1000 3/1000 0 -2/1000 3/1000 0 ...
% % 1 2/1000 3/1000 -1/1000 1.5/1000 0
% % 1 -2/1000 -3/1000 0 0/1000 -3/1000 0 0/1000 0/1000 0 -2/1000 0/1000 0 ...
% % 1 2/1000 3/1000 -1/1000 -1.5/1000 0
% % 1 0/1000 -3/1000 0 2/1000 -3/1000 0 2/1000 0/1000 0 0/1000 0/1000 0 ...
% % 1 2/1000 3/1000 1/1000 -1.5/1000 0];
num_elems = 30000;
x = 10;
y = num_elems/x;
rect = zeros(num_elems,19);
rect(:,1) = 1;
rect(:,14) = 1;
widthss = (2e-3 + 2e-3)/x;
lengthss = (3e-3 +3e-3)/y;
rect(:,15) = widthss;
rect(:,16) = lengthss;
val = 1;
for m = 0:x-1
    for n = 0:y-1
        rect(val,2:4) = [-2e-3 + m*widthss,-3e-3 + n*lengthss,0];
        rect(val,5:7) = [-2e-3 + widthss + m*widthss,-3e-3 + lengthss + n*lengthss,0];
        rect(val,8:10) = [-2e-3 + widthss + m*widthss,-3e-3 + n*lengthss,0];
        rect(val,11:13) = [-2e-3 + m*widthss,-3e-3 + lengthss + n*lengthss,0];
        rect(val,17:19) = [-2e-3 + m*widthss + widthss/2, -3e-3 + n*lengthss + lengthss/2, 0];
        val = val + 1;
    end
end
center=[0 0 0];
focus=[0 0 5e-3];
Th = xdc_rectangles (rect, center, focus);
ele_size = 1/1000;
xdc_show(Th)
points = [2.8e-3,3.75e-3,5e-3;2e-3,2.25e-3,5e-3;1.2e-3,2.25e-3,5e-3;1.2e-3,0,5e-3;0,0,5e-3];
x_a = [1.4,1,0.6,0.6,0];
y_b = [1.25,.75,.75,0,0];
[h, start_time] = calc_h(Th,points);
length_h = size(h,1);
t = start_time;
for m = 1:(length_h-1)
    t = [t; start_time+ele_size*m];
end
figure;
for k = 1:5
    subplot(5,1,k)
    plot(t,h(:,k));
    title(['Impulse Response at z/a = 10, b/a = 1.5, x/a = ', num2str(x_a(k)), ' and y/b = ',num2str(y_b(k))])
    xlabel('t')
    ylabel('h(x,t)')
    ylim([0 2000])
    hold on
end
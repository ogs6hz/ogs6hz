path(path,'/Users/Livvy/Documents/BME 8730/Field_II_ver_3_24_windows');
field_init(-1)
num_elems = 30000;
x = 10;
y = num_elems/x;
rect = zeros(num_elems,19);
rect(:,1) = 1;
rect(:,14) = 1;
widthss = (2e-3 + 2e-3)/x;
lengthss = (3e-3 +3e-3)/y;
R = ((2e-3)^2+(3e-3)^2);
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
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda = R/2.5; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/385)2.5 
f0 = c/lambda; % Transducer center frequency [Hz] 
func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
ampmag = [];
xpts = 0:0.1:1.5;
ypts = 0:0.1:1.5;
zpts = 0:0.1:1;
hp = [];
rvals = [];
for i = 1:length(xpts)
    fprintf('\n %d : ',i)
    for j = 1:length(ypts)
    for k = 1:length(zpts)
        fprintf(' %d ',k)                        % Short piece of code to provide record of progress
        if (k/20)==round(k/20)
            fprintf('\n')
        end
        [h, start_time_h(i,k)] = calc_h(Th,[i,j,k]);
        hconv = conv(h,func);
        hpconv = diff(hconv);
        ampmidpoint = round(size(hpconv,2)/2);
        ampmid = hpconv(ampmidpoint-round(ampmidpoint/8):ampmidpoint+round(ampmidpoint/8));
        r = round(sqrt(i^2+j^2));
        ampmag(r,k) = max(ampmid)-min(ampmid); 
        rvals(r) = r;
    end
    end
end
figure;
radial = rvals./R; %y = 0 at every point tested
axial = zpts./(R.^2)./lambda;
plot(radial,ampmag')
figure;
mesh(radial,axial,ampmag')
title('Figure 10')
xlabel('Radial Distance (units of a)')                        
ylabel('Axial Distance (units of a^2/lambda)')
zlabel('|p|')
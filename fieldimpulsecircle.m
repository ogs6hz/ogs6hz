path(path,'/Users/Livvy/Documents/BME 8730/Field_II_ver_3_24_windows');
%path(path,'/home/user/field_II/m_files');
field_init(-1)
%field_debug(1)
% Set initial parameters
R=1; % Radius of transducer
ele_size=1/1000; % Size of mathematical elements
% Define the transducer
Th = xdc_piston(R,ele_size);
points = [1.1987,.3,5;.954,.3,5;.8485,.3,5;.4,.3,5;.08,.03,5;0,0,5];
r_a = [1.2,1,0.9,0.5,0.1,0];
[h, start_time] = calc_h(Th,points);
length_h = size(h,1);
t = start_time;
for i = 1:(length_h-1)
    t = [t; start_time+ele_size*i];
end
figure;
for k = 1:6
%     [h, start_time] = calc_h(Th,points(k,:));
    subplot(6,1,k)
    plot(t,h(:,k));
    title(['Impulse Response at z/a = 5 and rho/a = ', num2str(r_a(k))])
    xlabel('t')
    ylabel('h(x,t)')
    ylim([0 2000])
    hold on
end
%% 
R=1; % Radius of transducer
ele_size=1/1000; % Size of mathematical elements
% Define the transducer
Th = xdc_piston(R,ele_size);
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda = R/2.5; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/385)2.5 
f0 = c/lambda; % Transducer center frequency [Hz] 
func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
ampmag = [];
xpts = 0:0.5:1.5;
ypts = 0:0.5:1.5;
zpts = 0:0.25:1;
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
title('Figure 7')
xlabel('Radial Distance (units of a)')                        
ylabel('Axial Distance (units of a^2/lambda)')
zlabel('|p|')
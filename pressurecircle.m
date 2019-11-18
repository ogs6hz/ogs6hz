clear all                                   % Clear work space from previous runs etc
zvalues = 490:5:510; %490:0.5:510;  
zsize = length(zvalues);
vel = 1540;                                   % Speed of sound - all units MKS
a = 100;
rvalues = 0:20:120; %0:5:120;
rsize = length(rvalues);
fs = 100e6; % Sampling frequency [Hz]
tvalues = 0.32:0.0001:0.36;
tsize = length(tvalues);
lambda = a/2.5; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/385)2.5 
f0 = vel/lambda; % Transducer center frequency [Hz] 
func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
% figure
ampmag = zeros(rsize,zsize);
for j = 1:zsize
    z = zvalues(j);
    fprintf('\n %d : ',j)
for i = 1:rsize
    fprintf(' %d ',i)                        % Short piece of code to provide record of progress
    if (i/20)==round(i/20)
            fprintf('\n')
    end
    r1 = rvalues(i);
    t1 = round(z/vel,4);
    t2 = round((1/vel)*sqrt(z^2 + (r1-a)^2),4);
    t3 = round((1/vel)*sqrt(z^2 + (r1+a)^2),4);
    h = [];
    for k = 1:tsize
        t = tvalues(k);
        if r1<a
            if t<t1 %define impulse response
                h(k) = 0*vel;
            elseif (t1<t)&&(t2>t)
                h(k) = 1*vel;
            elseif (t2<t)&&(t3>t)
                h(k) = (1/pi)*acos(((vel^2*t^2)-z^2+r1^2-a^2)/(2*r1*sqrt((vel^2*t^2)-z^2)))*vel;
            else
                h(k) = 0*vel;
            end
        else
            if t<t2
                h(k) = 0*vel;
            elseif (t2<t)&&(t3>t)
                h(k) = (1/pi)*acos(((vel^2*t^2)-z^2+r1^2-a^2)/(2*r1*sqrt((vel^2*t^2)-z^2)))*vel;
            else 
                h(k) = 0*vel;
            end
        end
    end
     velconv = conv(h,func); %convolve the sine wave with the impulse response
%     plot(velconv)
    press = diff(velconv); %take the derivative to get pressure
    figure
    plot(press)
    ampmidpoint = round(size(press,2)/2);
    ampmid = press(ampmidpoint-round(ampmidpoint/16):ampmidpoint+round(ampmidpoint/16));
    ampmag(i,j) = max(ampmid)-min(ampmid); %take the amplitude at the middle section of the wave (ignore the two ends)
end
end
figure;
radial = rvalues./a;
axial = zvalues./(a.^2)./lambda;
plot(radial,ampmag)
figure;
mesh(radial,axial,ampmag');
title('Figure 7')
xlabel('Radial Distance (units of a)')                        
ylabel('Axial Distance (units of a^2/lambda)')
zlabel('|p|')
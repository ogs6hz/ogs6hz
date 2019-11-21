clear all                                   % Clear work space from previous runs etc
% r0=0;                              % Number of array elements
% theta0=0;
% theta1 = 45;
% arc = sqrt(z^2 + r0^2 + r1^2 - 2*r0*r1*cos(theta1-theta0));
z = 500;  
vel = 1540;                                   % Speed of sound - all units MKS
a = 100;
rvalues = [120, 100, 90, 50, 10, 0];
rsize = length(rvalues);
tvalues = 0.32:0.00001:0.36;
tsize = length(tvalues);
figure
for i = 1:rsize
    r1 = rvalues(i);
    t1 = z/vel;
    t2 = (1/vel)*sqrt(z^2 + (r1-a)^2);
    t3 = (1/vel)*sqrt(z^2 + (r1+a)^2);
    h = [];
    for k = 1:tsize
        t = tvalues(k);
        if r1<a
            if t<t1
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
    subplot(rsize,1,i)
    plot(tvalues,h)
    title(['Impulse Response at z/a = 5 and rho/a = ', num2str(r1/a)])
    xlabel('t')
    ylabel('h(x,t)')
    ylim([0 2000])
    hold on
end
%% 
clear all
z = 15.4;  
vel = 1540;                                   % Speed of sound - all units MKS
avalues = 0:0.1:1.5;
lambdavalues=avalues./2.5; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/385)2.5 
asize = length(avalues);
r1 = 0;
tvalues = 0:1/100e6:0.07;
%tvalues = 0.0012; %t-value so t = t1
tsize = length(tvalues);
figure
amp = [];
ampmag = zeros(asize);
%syms x
for i = 1:asize
    fprintf('%d ',i)                        % Short piece of code to provide record of progress
    if (i/20)==round(i/10)
            fprintf('\n')
    end
    a = avalues(i);
    lambda = lambdavalues(i);
    t1 = z/vel;
    t2 = round((1/vel)*sqrt(z^2 + (r1-a)^2),5);
    t3 = round((1/vel)*sqrt(z^2 + (r1+a)^2),5);
    dh = [];
    func = [];
    for k = 1:tsize
        t = tvalues(k);
        if r1<a
            if t<t1
                dh(k) = 0*vel;
            elseif t == t1
                dh(k) = 1*vel; %derviative at t1 = change in h over change in time = (1*vel-0*vel)/1 = 1*vel
            elseif (t1<t)&&(t2>t)
                dh(k) = 0*vel;
            %elseif (t2<t)&&(t3>t)
                 %h(k) = (1/pi)*acos(((vel^2*t^2)-z^2+r1^2-a^2)/(2*r1*sqrt((vel^2*t^2)-z^2)))*vel;
            elseif t == t2 %for on-axis, t2 = t3 and change at t2 = (0*vel-1*vel)/1 = -1*vel
                dh(k) = -1*vel;
            elseif t>t2
                dh(k) = 0*vel;
            end
        else
            if t<t2
                dh(k) = 0*vel;
            %elseif (t2<t)&&(t3>t)
                %dh(k) = (1/pi)*acos(((vel^2*t^2)-z^2+r1^2-a^2)/(2*r1*sqrt((vel^2*t^2)-z^2)))*vel;
            else 
                dh(k) = 0*vel;
            end
        end
    end
    %rho = 1060; %density of tissue
    f0=vel/lambda; % Transducer center frequency [Hz] 
    fs=100e6; % Sampling frequency [Hz]
    func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
    amp = conv(dh,func);
    amp = abs(amp); 
    ampmag(i) = max(amp);
    %plot(func)
    %gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
end
figure;
plot(avalues,ampmag)
axial = (avalues.^2)./lambdavalues;
%mesh(avalues,axial,ampmag);
title('Figure 7')
xlabel('Radial Distance (units of a)')                        % Be sure to account for mm or m here
%ylabel('Axial Distance (units of a^2/lambda)')
%zlabel('|p|')
    %% 
    
clear all
z = 15.4; 
vel = 1540;                                   % Speed of sound - all units MKS
a = 10;
lambda=4; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/385)2.5 
r1 = 0;
tvalues = 0.00199:0.000001:0.01299;
%tvalues = 0.0012; %t-value so t = t1
tsize = length(tvalues);
figure
%syms x
    t1 = z/vel;
    t2 = round((1/vel)*sqrt(z^2 + (r1-a)^2),5);
    t3 = round((1/vel)*sqrt(z^2 + (r1+a)^2),5);
    for k = 1:tsize
        t = tvalues(k);
        if r1<a
            if t<t1
                dh(k) = 0*vel;
            elseif t == t1
                dh(k) = 1*vel; %derviative at t1 = change in h over change in time = (1*vel-0*vel)/1 = 1*vel
            elseif (t1<t)&&(t2>t)
                dh(k) = 0*vel;
            %elseif (t2<t)&&(t3>t)
                 %h(k) = (1/pi)*acos(((vel^2*t^2)-z^2+r1^2-a^2)/(2*r1*sqrt((vel^2*t^2)-z^2)))*vel;
            elseif t == t3 %for on-axis, t2 = t3 and change at t2 = (0*vel-1*vel)/1 = -1*vel
                dh(k) = -1*vel;
            elseif t>t3
                dh(k) = 0*vel;
            end
        else
            if t<t2
                dh(k) = 0*vel;
            %elseif (t2<t)&&(t3>t)
                %dh(k) = (1/pi)*acos(((vel^2*t^2)-z^2+r1^2-a^2)/(2*r1*sqrt((vel^2*t^2)-z^2)))*vel;
            else 
                dh(k) = 0*vel;
            end
        end
    end
    %rho = 1060; %density of tissue
    f0=vel/lambda; % Transducer center frequency [Hz] 
    fs=100e6; % Sampling frequency [Hz]
    func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
    plot(func)
    %gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
    amp = conv(dh,func);
    amp = abs(amp);
    amp = amp/max(amp);%normalize pressure
    ampsize = size(amp,2);
    xaxis = 0:a/(ampsize-1):a;
    xaxis = xaxis/max(xaxis);
    figure
    plot(dh)
    figure
    plot(xaxis,amp)
%% 
clear all
z = 154;  
vel = 1540;                                   % Speed of sound - all units MKS
a = 70.5716657;
r1 = 0.0016;
tvalues = 0.0999:0.0001:0.1111;
tsize = length(tvalues);
figure
t1 = z/vel;
t2 = (1/vel)*sqrt(z^2 + (r1-a)^2);
t3 = (1/vel)*sqrt(z^2 + (r1+a)^2);
dh = [];
rho = 1060; %density of tissue
for i = 1:tsize
    t = tvalues(i);
    d1 = dirac(t-(z/vel));
    d2 = dirac(t-(sqrt(z^2+a^2))/vel);
    dh(i) = rho*vel*(d1-d2);
end
%trans = conv(2*rho*vel*exp((-sqrt(z^2+a^2)*s)/vel),sinh((sqrt(z^2+a^2)-z)/vel*s);
f0=3e6; % Transducer center frequency [Hz] 
lambda=vel/f0; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/3e6)2.5 
fs=100e6; % Sampling frequency [Hz]
func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
plot(func)
%gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
amp = conv(dh,func);
figure
plot(dh)
figure
plot(amp)
%% 
z = 10000; %z distance 
b = 150; %radius of length
vel = 1540;  %speed of sound
a = 100; %radius of width
xvalues = [140, 100, 60, 60, 0]; %x-direction of field point
yvalues = [187.5, 112.5, 112.5, 0, 0]; %y-direction of field point
xsize = length(xvalues);
tvalues = -0.001:0.0000001:0.005; %-0.001:0.00001:0.005; %6.4934:0.0001:6.5;
tsize = length(tvalues);
figure
for i = 1:xsize
    x = abs(xvalues(i));
    y = abs(yvalues(i));
    h = [];
    arc = sqrt(z^2+x^2+y^2);
    if abs(x) <= a && abs(y) <= b
        small = [a+x,a-x,a+x,a-x];
        elmall = [b+y,b+y,b-y,b-y];
        coeffall = [1,1,1,1];
    elseif abs(x) > a && abs(y) > b
        small = [2*a,2*a,x-a,x-a];
        elmall = [2*b,y-b,2*b,y-b ];
        coeffall = [1,-1,-1,-1];
    elseif abs(x) > a && abs(y) <= b 
        small = [2*a,2*a,x-a,x-a];
        elmall = [b+y,b-y,b+y,b-y];
        coeffall = [1,1,-1,-1];
    elseif abs(x) <= a && abs(y) > b 
        small = [a+x,a+x,a-x,a-x];
        elmall = [2*b,y-b,2*b,y-b];
        coeffall = [1,-1,1,-1];
    end
    for k = 1:tsize
        tau = tvalues(k);
        t = tau + arc/vel; 
        t1 = z/vel;
        summation = 0;
        for m = 1:4
            sm = small(m);
            elm = elmall(m);
            coeff = coeffall(m);
            t2 = (1/vel)*sqrt(z^2 + sm^2);
            t3 = (1/vel)*sqrt(z^2 + elm^2);
            t4 = (1/vel)*sqrt(z^2 + sm^2 + elm^2);
            if t<t1
                var = 0;
            elseif (t1<=t)&&(t2>t)
                var = pi/2;
            elseif (t2<=t)&&(t3>t)
                var = pi/2 - acos(sm/sqrt((vel.^2)*(t.^2)-(z^2)));
            elseif (t3<=t)&&(t4>=t)
                var = pi/2 - acos(sm/sqrt((vel.^2)*(t.^2)-(z^2))) - acos(elm./sqrt((vel.^2).*(t.^2)-(z.^2)));
            else
                var = 0;
            end                
%             t14 = integral(@(t) dirac(t-arc/vel)*(pi/2), t1, t4);
%             t24 = integral(@(t) dirac(t-arc/vel).*acos(sm./sqrt((vel.^2).*((t-(t-arc/vel)).^2)-(z.^2))),t2,t4);
%             t34 = integral(@(t) dirac(t-arc/vel).*acos(elm./sqrt((vel.^2).*((t-(t-arc/vel)).^2)-(z.^2))),t3,t4);
            if var ~= 0
                summation = summation + coeff*var; %(t14 - t24 - t34);
            end
        end
        h(k) = (vel/(2*pi))*summation;  
    end
    subplot(xsize,1,i)
    plot(tvalues+arc/vel,h)
    title(['Impulse Response at z/a = 10, x/a = ', num2str(x/a), ' and y/b = ', num2str(y/b)])
    xlabel('t')
    ylabel('h(x,t)')
    %axis tight
    xlim([6.493 6.5])
    ylim([0 1540])
    hold on
end

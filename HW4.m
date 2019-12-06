z_foc=50e-3;                                % Range direction focal distance
num_elems=128;                              % Number of array elements
fs=100e6; % Sampling frequency [Hz]
pitch = 0.075e-3;
vel=1540; % Speed of sound [m/s]
f = 100e6;
lambda = vel/f; %lambda = v/f want a/lambda = 2.5, so a/(v/f)=2.5 a = (1540/385)2.5 
f0 = vel/lambda; % Transducer center frequency [Hz] 
func = sin(2*pi*f0*(0:1/fs:10/f0)); %sin wave of 10 cycles
w=2*pi*f0*(0:1/fs:10/f0);                                   % Angular frequency radians
x_pts=[-10:1:10].*1e-3;                     % Define X-direction field locations
z_pts=[25:2:75].*1e-3;               % Define Z-direction field locations
tstep=1./max(f0*(0:1/fs:10/f0));                            % Time steps after using Inverst FFT
ns=length(f0*(0:1/fs:10/f0));                               % Number of samples
t=[1:ns].*tstep;                            % Define time axis

for m=1:length(t)                       % Loop over field locations
    fprintf('%d ',j)                        % Short piece of code to provide record of progress
    if (j/20)==round(j/20)
            fprintf('\n')
    end
    for i=1:length(x_pts) 
        sum_wave=0;  % Initialize sum of waveforms to zero before starting loop
        for k=1:num_elems
            x_elem(k)=((k-1)-(num_elems-1)./2).*pitch;  % Calculate locations of each array element
            tdel = sqrt((x_elem(k)-x_pts(i))^2)/vel; %apply time delay
            wave=func(t(m)).*exp(-j*w*tdel);    % Apply time delay so 0 for t<0
            sum_wave = sum_wave+wave; 
            %sum_wave should be a matrix with a different sum for each
            %point
        end
        sum_pulse_t=ifft(sum_wave);                  % Convert to time domain
        field_val(m,i) = max(sum_pulse_t)-min(sum_pulse_t);
        %field_val(m,i)=max(abs(hilbert(sum_pulse_t)));     % Field values in beamplot - Hilbert finds envelope of sinsoidal fn
    end
    field_db(m,:)=real(field_val(m,:)./max(field_val(m,:)));     % Normalize by peak value at each range
    %field_db(m,:)=20*log10(field_val(m,:));                 % Convert to dB
end

figure
mesh(x_pts*1e3,z_pts*1e3,field_db)
%axis([-2 2 -30 0])
title('Beamplot at z = 50 mm')                          % Label according to question
xlabel('Range (mm)')                        % Be sure to account for mm or m here
ylabel('Azimuth (mm)')
zlabel('Sum of Sine Waves'); 

figure
scats = 0 + (10-0)*rand(length(x_pts),length(z_pts)); %random scatters from 0 to 10
psfconv = conv2(field_db,scats);
imagesc(psfconv)
title('Convoluted Func at z = 50 mm')                          % Label according to question
xlabel('Range (mm)')                        % Be sure to account for mm or m here
ylabel('Azimuth (mm)')
zlabel('Sum of Sine Waves'); 
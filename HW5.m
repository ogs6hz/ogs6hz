close all
sig33 = 1200 * 8.8543e-12;
rho  = 7820;
cE11 = 1.3745e11; 
cE33 = 1.2574e11;                 
cE12 = 0.8790e11;                 
cE13 = 0.9230e11;                 
e31 = -9.440;                   
e33 =  22.495; 
thick = 0.012; %elevation
width= 0.0003; %azimuth
len = 0.0008;  %thickness
area = width*thick;
Z1 = 1.5e6*area;
%Z2 = 4e6*area;
Z2 = 10e6*area;

cE33 = cE33*(1-cE13^2)/(cE11*cE33); %modified piezoelectric coeff
e33 = e33-(e31*cE13)/cE11; %modified elastic stiffness

fs=1e2;
fmax=40e6;
ns=2^(nextpow2(fmax/fs));
fmax=ns*fs;
f=(1:ns)*fs;
ts=1/fmax;
t=(1:ns)*ts;
j=sqrt(-1);
% 2 cycle sq wave. This makes the fundamental frequency more evident.
h1 = heaviside(t)-heaviside(t-2.5e-7); %half cycle
%h2 = heaviside(t-2.5e-7)-heaviside(t-5e-7);
h2 = heaviside(t-5e-7)-heaviside(t-7.5e-7);
h=h1+h2; 
fsq = (fft(h));
w=2*pi*f;                                   % Angular frequency radians
ns = length(fsq);

h = e33/sig33;
cD = cE33*(1+e33^2/(cE33*sig33));
Z0 = sqrt(rho*cD);
B = w.*sqrt(rho/cD);
c0 = sig33*area/len;
Zc = Z0*area;
%% 
%Transmitter
A = zeros(3,3,ns);
A(1,1,:) = -j.*((Zc.*cot(B.*len))+(Z1./-j));
A(1,2,:) = -j.*(Zc.*csc(B.*len));
A(1,3,:) = -j.*(h./w);
A(2,1,:) = -j.*(Zc.*csc(B.*len));
A(2,2,:) = -j.*((Zc.*cot(B.*len))+(Z2./-j));
A(2,3,:) = -j.*(h./w);
A(3,1,:) = -j.*(h./w);
A(3,2,:) = -j.*(h./w);
A(3,3,:) = -j.*(1./(w.*c0));

V3 = Zel/(Zel+50)
v3 = 1;
b = [0;0;v3];
 for i = 1:ns
x(:,i) = A(:,:,i)\b; %Solve for V1 & V2 (V3 = 1)
 end
v1 = x(1,:).*w;
v2 = x(2,:).*w;
F1 = -Z1*v1; %F1 = -z1v1 : z1 is front material
F1sq = F1.*fsq;
F2 = -Z2*v2; %F2 = -z2v2 : z2 is rear material
I3 = x(3,:); %Move back to F1 and F2, and solve for I3
Zel = real(I3/v3); %Ze = V3/I3
Zelsq = Zel.*fsq;
F1_t= ifft(F1);
F1_t = F1_t./max(F1_t);
F1sq_t = ifft(F1sq);
F1sq_t = F1sq_t./max(F1sq_t);
Zel_t = ifft(Zel);
Zelsq_t = ifft(Zelsq);
F1 = real(20*log10(abs(F1./F1(fmax/fs))));
F1sq = real(20*log10(abs(F1sq./F1sq(fmax/fs))));
% 
% bw = .8;
% sig=bw*fc/100;                              % Width of Gaussian
% gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
% tdel = 1;
% gauss_pulse=gauss_pulse.*exp(-j*w*tdel);    % Apply time delay so 0 for t<0

figure
plot(f,F1)
xlabel('Frequency (Hz)')
ylabel('F1')

figure
plot(f,Zel)
xlabel('Frequency (Hz)')
ylabel('|Zel|')


% tstep=1./max(f);                            % Time steps after using Inverst FFT
% t=[1:ns].*tstep;                            % Define time axis 

figure
plot(real(F1_t(1:500)))
xlabel('Time (s)')
ylabel('F1')

figure                                    % Plot waveform as a diagnostic
plot(real(Zel_t(1:500)))
xlabel('Time (s)')
ylabel('|Zel|')

figure
plot(f,F1sq)
xlabel('Frequency (Hz)')
ylabel('F1 2 Cycle Square Wave')

figure
plot(real(F1sq_t(1:500)))
xlabel('Time (s)')
ylabel('F1 2 Cycle Square Wave')
%ifft(multiply Fs by gaussian/square wave) to get time fft 
%Plot this Ze vs frequency
%Plot phase vs frequncy

figure
plot(f,Zelsq);
xlabel('Frequency (Hz)')
ylabel('|Zel| 2 Cycle Square Wave');

figure
plot(real(Zelsq_t(1:500)));
xlabel('Time(s)')
ylabel('|Zel| 2 Cycle Square Wave');
%% 
cE33mod = cE33*(1-cE13^2)/(cE11*cE33); %modified piezoelectric coeff
e33mod = e33-(e31*cE13)/cE11; %modified elastic stiffness

%% 
%Receiver
A = zeros(3,3,ns);
A(1,1,:) = 0;
A(1,2,:) = -j*(Zc/sin(B*len));
A(1,3,:) = -j*(h/w);
A(2,1,:) = 0;
A(2,2,:) = -j*((Zc/tan(B*len))+(Zc/-j));
A(2,3,:) = -j*(h/w);
A(3,1,:) = -v3;
A(3,2,:) = -j*(h/w);
A(3,3,:) = -j*(1/(w*c0));

v1 = -F1/z1;
b = [j*(Zc/tan(B*len))+Zc;j*(Zc/sin(B*len));j*(h/w)];
for i = 1:size(w,1)
    x(:,i) = A(:,:,i)\b; %Solve for V2 & V3 (V1 = 1)
end
v3 = x(1,:);
v2 = x(2,:);
I3 = x(3,:); %Move back to F1 and F2, and solve for I3
Zel = I3/v3; %Ze = V3/I3

bw = 30;
sig=bw*fc/100;                              % Width of Gaussian
gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
tdel = 1;
gauss_pulse=gauss_pulse.*exp(-j*w*tdel);    % Apply time delay so 0 for t<0

ff_t = ifft(Zel.*gauss_pulse);

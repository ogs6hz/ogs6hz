sig0 = 8.85e-12;
cE11 = 12.1e10;
cE12 = 7.54e10;
cE13 = 7.52e10;
cE33 = 11.1e10;
e31 = -5.4;
e33 = 15.8;
sig33sig0 = 830;
sig33 = sig33sig0*sig0;
rhoc = 7750;
k33 = 70;
kt = 49;
c11 = 0.53e10;
c12 = 0.31e10;
rhop = 1100;

v = 0:100;
cE33 = v.*(cE33-(2.*v.*(cE13-cE12).^2)/(v.*(cE11+cE12)+v.*(cE11+cE12)))+v.*cE11;
e33 = v.*(e33+(2.*v.*e31.*(cE13-cE12))/(v.*(cE11+cE12)+v.*(cE11+cE12)));
sig33 = (v.*(sig33+(2.*v.*(e31)^2)/(v.*(cE11+cE12)+v.*(cE11+cE12))))*sig0;
sig33sig0 = sig33/sig0;
cD33 = cE33+(e33).^2./sig33;
h33 = e33./sig33;
B33 = 1./sig33;
rho = v.*rhoc+v.*rhop;
kt = e33./sqrt(cD33.*sig33);
Z = sqrt(cD33.*rho);
vel = sqrt(cD33./rho);

figure
subplot(4,2,1);
plot(v,rho);
xlabel('Volume Fraction Ceramic (%)')
ylabel(' rho (g/cm^3)')
subplot(4,2,2)
plot(v,kt);
xlabel('Volume Fraction Ceramic (%)')
ylabel(' Kt (%)')
subplot(4,2,3)
plot(v,sig33sig0);
xlabel('Volume Fraction Ceramic (%)')
ylabel('Epsilon33/Epsilon0')
subplot(4,2,4)
plot(v,vel);
xlabel('Volume Fraction Ceramic (%)')
ylabel('Vl (m/sec)')
subplot(4,2,5)
plot(v,cD33);
xlabel('Volume Fraction Ceramic (%)')
ylabel('c^D33 (10^10 N/m^2)')
subplot(4,2,6)
plot(v,Z);
xlabel('Volume Fraction Ceramic (%)')
ylabel('Z (Mrayls)')
subplot(4,2,7)
plot(v,e33);
xlabel('Volume Fraction Ceramic (%)')
ylabel('e33 (C/m^2)')
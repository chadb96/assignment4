close all;
clear;

%constants
R1=1;
G1=1/R1;
c=0.25;
R2=2;
G2=1/R2;
L=0.2;
R3=10;
G3=1/R3;
a=100;
R4=0.1;
G4=1/R4;
Ro=1000;
Go=1/Ro;

G= [G1 0 0 0 0 0 0;
    -G1 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -a 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+Go];
C= [0 0 0 0 0 0 0;
    -c c 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];
V= zeros(7,1);
F= zeros(7,1);
f1= figure;
f2= figure;
f3= figure;
f4= figure;

for v=-10:0.5:10
    F(1,1)= v;
    V= G\F;
    set(0, 'CurrentFigure', f1)
    scatter(v, V(7,1), 'r.')
    hold on
    title('Vo in DC')
    
    set(0, 'CurrentFigure', f2)
    scatter(v, V(4,1), 'b.')
    hold on
    title('V3 in DC')
end

w=logspace(1,2,500);
F(1)=1;
for i=1:length(w)
    Vac = (G+C*i*i)\F;          
    set(0, 'CurrentFigure', f2)
    semilogx(w(i), abs(Vac(7,1)), 'g.')
    hold on
    title('Vo in AC case')
    
    dB = 20*log(abs(Vac(7,1))/F(1));   
    set(0, 'CurrentFigure', f3)
    plot(i, dB, 'c.')
    hold on
    title('Gain with respect to Frequency')
    xlabel 'Frequency (Hz)';
    ylabel 'Vout';
end

% Calculating and plotting AC sweep 
cs = 0.25+0.05.*randn(1,1000);
Vgain = zeros(1000,1);
for j = 1:length(Vgain)
    c = cs(j);
    C(2,1) = -c;
    C(2,2) = c;
    Vac = (G+C*1j*pi)\F;                 
    Vgain(j,1) = abs(Vac(7,1))/F(1);    
end

% Histogram 
set(0, 'CurrentFigure', f4)
hist(Vgain,50);
title 'Histogram of Voltage Gain';

%Part 2
step= 0.001;
A= C/step + G; 
Vold= [1;0;0;0;0;0;0];
Vnew= [1;0;0;0;0;0;0];
count= 1;
for i=0.001:0.001:1000*10^-3
    if count<30
            F(1,1)= 0;
            V= G\F;
            Vout(count,1)= i;
            Vout(count,2)= V(7);    
    else
            F(1,1)= 1;
            Vnew= inv(A)*(F+C*Vold./step);
            Vold= Vnew;
            Vout(count,1)= i;
            Vout(count,2)= Vnew(7);
            Vout(count,3)= Vnew(1);   
    end    
            count= count + 1;
end
% Step Function

figure (5);
plot(Vout(:,1),Vout(:,2),'r');
hold on;
plot(Vout(:,1),Vout(:,3),'b');
title 'Step Function Plot';
xlabel 'Time (s)';
ylabel 'Vout';
legend ("Output","Input");
grid on;

% Frequency Plot
VoutFnoshift= fft(Vout(:,2));
VinFnoshift= fft(Vout(:,3));
VinF= fftshift(VinFnoshift);
VoutF= fftshift(VoutFnoshift);

figure(6);
semilogy (abs(VoutF));
hold on;
semilogy(abs(VinF));
legend("Output","Input");
title 'Step Function in Frequency';
xlabel 'Frequency (Hz)';
ylabel 'Magnitude';
grid on;

count= 1;
step= 0.001;
Vold= [0;0;0;0;0;0;0];
Vnew= [1;0;0;0;0;0;0];
%Part 2 B
for i=0.001:0.001:1000*10^-3
            F(1,1)= sin(2*pi*(1/0.03)*i);
            Vnew= inv(A)*(F + C*Vold./step);
            Vold= Vnew;
            Voutb(count,1)= i;
            Voutb(count,2)= Vnew(7);
            Voutb(count,3)= F(1,1);
            count= count+1;
end
figure (7);
plot(Voutb(:,1), Voutb(:,2),'r');
hold on;
plot(Voutb(:,1), Voutb(:,3),'b');
title 'Sine Function Input';
xlabel 'Time (s)';
ylabel 'Vout';
legend ("Output","Input");
grid on;

% Frequency Plot
VoutFbnoshift= fft(Voutb(:,2));
VinFbnoshift= fft(Voutb(:,3));
VoutFb= fftshift(VoutFbnoshift);
VinFb= fftshift(VinFbnoshift);

figure(8);
semilogy (abs(VoutFb));
hold on;
semilogy(abs(VinFb));
legend("Output","Input");
title 'Sine Function in Freq Domain';
xlabel 'Frequency (Hz)';
ylabel 'Magnitude';
grid on;

%Part 2 C
count= 1;
step= 0.001;
Vold= [0;0;0;0;0;0;0];
Vnew= [1;0;0;0;0;0;0];

for i=0.001:0.001:1000*10^-3

                F(1,1)= (exp(-(i-0.06).^2/(2*step)));
                Vnew= inv(A)*(F+C*Vold./step);
                Vold= Vnew;
                Voutc(count,1)= i;
                Voutc(count,2)= Vnew(7);
                Voutc(count,3)= F(1,1);
                count= count+1;
end
figure (9);
plot(Voutc(:,1), Voutc(:,2),'r');
hold on;
plot(Voutc(:,1), Voutc(:,3),'b');
title 'Gaussian Pulse Input';
xlabel 'Time (s)';
ylabel 'Vout';
legend ("Output","Input");
grid on;

% Frequency Plot
VoutFcnoshift= fft(Voutc(:,2));
VinFcnoshift= fft(Voutc(:,3));
VinFc= fftshift(VinFcnoshift);
VoutFc= fftshift(VoutFcnoshift);

figure(10);
semilogy (abs(VoutFc));
hold on;
semilogy(abs(VinFc));
legend("Output","Input");
title 'Gaussian Pulse in Frequency';
xlabel 'Frequency (Hz)';
ylabel 'Magnitude';
grid on;


Cn= 0.00001;
C= zeros(7,7);
C(2,1)= -c;
C(2,2)= c;
C(3,3)= -L;
C(4,4)= -Cn;
C(6,4)= -Cn;
F = zeros(7,1);

%Part 3
count= 1;
step= 0.01;
Vold= [1;0;0;0;0;0;0];
Vnew= [1;0;0;0;0;0;0];

for i=0.001:0.001:1000*10^-3
                F(7,1)= randn(1,1);
                F(1,1)= (exp(-(i-0.06).^2/(2*step)));
                Vnew= inv(A)*(F + C*Vold./step);
                Vold= Vnew;
                Vout3(count,1)= i;
                Vout3(count,2)= Vnew(7);
                Vout3(count,3)= F(1,1); 
                count= count+1;
end
% Vout with noise
figure (11);
plot(Vout3(:,1), Vout3(:,2),'r');
hold on;
plot(Vout3(:,1), Vout3(:,3),'b');
title 'Vout with Input In as Noise';
xlabel 'Time (s)';
ylabel 'Vout';
legend ("Output","Input");
grid on;

% Frequency Plot
Vout3noshift= fft(Vout3(:,2));
Vin3noshift= fft(Vout3(:,3));
VinF3= fftshift(Vin3noshift);
VoutF3= fftshift(Vout3noshift);

figure(12);
semilogy (abs(VoutF3));
hold on;
semilogy(abs(VinF3));
legend("Output","Input");
title 'Vout in Frequency with Noise';
xlabel 'Frequency (Hz)';
ylabel 'Magnitude';
grid on;

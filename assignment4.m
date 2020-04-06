function [] = assignment4()
clear
clc
close all
warning off

addpath('C:\Users\eshlden\Documents\Personal\School\PHYS 4007 - 4th Year Lab\MatlabLaboratoryHelp');

myCircuit = initialize();
q1(myCircuit);
q2(myCircuit);
q3(myCircuit);
end

function [myCircuit] = initialize()
myCircuit = circ();
myCircuit = myCircuit.add('vin',0,1,1);
myCircuit = myCircuit.add('r1',1,2,1);
myCircuit = myCircuit.add('c',1,2,0.25);
myCircuit = myCircuit.add('r2',2,0,2);
myCircuit = myCircuit.add('l',2,3,0.2);
myCircuit = myCircuit.add('r3',3,0,determineR3('IV.mat'));
myCircuit = myCircuit.add('b',0,4,100,'r3');
myCircuit = myCircuit.add('r4',4,5,0.1);
myCircuit = myCircuit.add('ro',5,0,1000);
myCircuit = myCircuit.setIn('vin');
myCircuit = myCircuit.setup;
end

function [] = q1(myCircuit)
myCircuit.printMat('G',false);
myCircuit.printMat('C',false);
myCircuit.printEqn();


myCircuit = myCircuit.DCsweep(linspace(-10,10,100));
f = myCircuit.output(5,3);
createFigure(f,[pwd,'\'],'DC Sweep of Input Voltage','q1dc');

myCircuit = myCircuit.ACsweep(linspace(0,100,101));
f = myCircuit.output(5);
createFigure(f,[pwd,'\'],'AC Sweep Frequency at Input Voltage 1V','q1ac');

f = varyC(myCircuit);
createFigure(f,[pwd,'\'],'Impacts of Gaussian Variation of Capacitanc','q1c');
end

function [f] = varyC(myCircuit)
trials = 1000;
cValues = zeros(1,trials);
f = 0.5; 
Vin = 1;
gain = zeros(1,trials);

for n = 1:trials
    cValues(n) = normrnd(0.25,0.05);
    myCircuit = myCircuit.change('c',cValues(n));
    myCircuit = myCircuit.setup();
    myCircuit = myCircuit.opPoint(Vin,f);
    gain(n) = 20*log10(abs(myCircuit.getSol(5)./Vin));
end

f = figure;
subplot(2,1,1)
plot(cValues,gain,'.');
xlabel('Capacitance (F)');
ylabel('Gain (dB)');
grid on;
title('Maximum Circuit Gain as a Function of Capacitance');

subplot(2,1,2)
histogram(gain);
xlabel('Gain (dB)');
ylabel('Count');
grid on;
title('Distribution of Gain with Variation of Capacitance');
end

function [] = q2(myCircuit)
t = linspace(0,1,1024);
V = zeros(1,length(t));
V(t>=0.03) = 1;
myCircuit = myCircuit.transient(t,V);
h = myCircuit.output(1,5);
createFigure(h,[pwd,'\'],'Step Input Transient Simulation','q2a');

f = 1/0.03;
V = sin(2*pi*f*t);
myCircuit = myCircuit.transient(t,V);
h = myCircuit.output(1,5);
createFigure(h,[pwd,'\'],'Sinusoidial Input Transient Simulation','q2b');

V = exp(-0.5*((t - 0.06)/0.03).^2);
myCircuit = myCircuit.transient(t,V);
h = myCircuit.output(1,5);
createFigure(h,[pwd,'\'],'Gaussian Input Pulse Transient Simulation','q2c');
end

function [] = q3(myCircuit)
Imean = 0;
Istd = 0.001;
myCircuit = myCircuit.add('in',3,0,Imean);
myCircuit = myCircuit.add('cn',3,0,0.00001);
myCircuit = myCircuit.setup;
myCircuit.printMat('C',true);

t = linspace(0,1,32768);
V = exp(-0.5*((t - 0.06)/0.03).^2);
myCircuit = myCircuit.transient(t,V,'in',Imean,Istd);
h = myCircuit.output(1,5);
createFigure(h,[pwd,'\'],'Gaussian Input Pulse Transient Simulation with noise source','q3a');

c = [0.0001,0.001,0.01];
for n = 1:length(c)
    t = linspace(0,1,1024);
    V = exp(-0.5*((t - 0.06)/0.03).^2);
    myCircuit = myCircuit.change('cn',c(n));
    myCircuit = myCircuit.transient(t,V,'in',Imean,Istd);
    h = myCircuit.output(1,5);
    subplot(2,1,1)
    title(sprintf('Transient Sweep (cn = %.2f mF)',c(n)*1000));
    subplot(2,1,2)
    title(sprintf('Frequency Spectrum (cn = %.2f mF)',c(n)*1000));
    createFigure(h,[pwd,'\'],sprintf('Gaussian Input Pulse Transient Simulation (cn = %.2f mF)',c(n)*1000),['q3b',num2str(n)]);
end
myCircuit = myCircuit.change('cn',0.00001);
N = [1024,32];
for n = 1:length(N)
    t = linspace(0,1,N(n));
    V = exp(-0.5*((t - 0.06)/0.03).^2);
    myCircuit = myCircuit.transient(t,V,'in',Imean,Istd);
    h = myCircuit.output(1,5);
    dt = 1000*(max(t) - min(t))/length(t);
    subplot(2,1,1)
    title(sprintf('Transient Sweep (Time step = %.2f ms)',dt));
    subplot(2,1,2)
    title(sprintf('Frequency Spectrum (Time step = %.2f ms)',dt));
    createFigure(h,[pwd,'\'],sprintf('Gaussian Input Pulse Transient Simulation (Time step = %.2f ms)',dt),['q3c',num2str(n)]);
end
end


function [R3] = determineR3(file,createPlot)
    if nargin == 1
       createPlot = false; 
    end
    
    load(file,'I','V');
    fitobject = fit(V',I',fittype('a*x'));
    R3 = 1./(fitobject.a);
    f = @(V) V./R3;
    
    fprintf('The resistance of R3 is %f Ohm\n\n',R3);
    
    if ~createPlot
        return;
    end
    
    figure;
    plot(V,I,'.-');
    hold on
    fplot(f,[min(V),max(V)]);
    grid on;
    ylabel('Current (A)');
    xlabel('Voltage (V)');
    title('IV Curve of Moddled Device');
    ylim([0,max(I)]);
end
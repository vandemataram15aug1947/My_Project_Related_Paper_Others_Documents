% PLL Modelling

Fs       =  50000;          % Sampling frequency = 50Khz
GridFreq =  50;             % Nominal Grid Frequency in Hz
Tfinal   =  0.2;            % Time the simulation is run for = 0.5 seconds
Ts       =  1/Fs;           % Sampling Time = 1/Fs
t        =  0:Ts:Tfinal;    % Simulation Time vector
wn       =  2*pi*GridFreq;  % Nominal Grid Frequency in radians

u  = zeros(size(t));
u1 = zeros(size(t));

% CASE 1 : Phase Jump at the Mid Point

 L=length(t);
 for n=1:floor(L)
 u(n)=sin(2*pi*GridFreq*Ts*n);
 end
 for n=1:floor(L)
 u1(n)=sin(2*pi*GridFreq*Ts*n);
 end
 for n=floor(L/2):L
 u(n)=sin(2*pi*GridFreq*Ts*n+pi/2);
 end

% CASE 2 : Harmonics

% L=length(t);
% for n=1:floor(L)
% u(n)=0.9*sin(2*pi*GridFreq*Ts*n)+0.1*sin(2*pi*5*GridFreq*Ts*n);
% end
% for n=1:floor(L)
% u1(n)=sin(2*pi*GridFreq*Ts*n);
% end

% CASE 3 : Frequency Shift

% L=length(t);
% for n=1:floor(L)
% u(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=1:floor(L)
% u1(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=floor(L/2):L
% u(n)=sin(2*pi*GridFreq*1.1*Ts*n);
% end

% CASE 4: Amplitude Variations

% L=length(t);
% for n=1:floor(L)
% u(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=1:floor(L)
% u1(n)=sin(2*pi*GridFreq*Ts*n);
% end
% for n=floor(L/2):L
% u(n)=0.8*sin(2*pi*GridFreq*Ts*n);
% end;

%Declare arrays used by the PLL process
Upd         =  zeros(size(t));
ynotch      =  zeros(size(t));
ynotch_buff =  zeros(size(t));
ylf         =  zeros(size(t));
SinGen      =  zeros(size(t));
Plot_Var    =  zeros(size(t));
Mysin       =  zeros(size(t));
Mycos       =  zeros(size(t));
theta       =  zeros(size(t));
werror      =  zeros(size(t));

% Notch filter design
c1 = 0.1;
c2 = 0.00001;
X  = 2*c2*wn*2*Ts;
Y  = 2*c1*wn*2*Ts;
Z  = wn*2*wn*2*Ts*Ts;
B_notch = [1 (X-2) (-X+Z+1)];
A_notch = [1 (Y-2) (-Y+Z+1)];

Mycos(2) = 1.0;
Mycos(1) = 1.0;

% Simulate the PLL process
for n=2:Tfinal/Ts 
    
    % Phase Detect
    Upd(1)    =  u(n)*Mycos(2);
    
    % Notch Filter
    ynotch(1) = -A_notch(2)*ynotch(2)-A_notch(3)*ynotch(3)+B_notch(1)*Upd(1)+B_notch(2)*Upd(2)+B_notch(3)*Upd(3);
   
    % Update the Upd array for future sample
    Upd(3)=Upd(2);
    Upd(2)=Upd(1);
    
    %---------------------------------------------------------------------%
    % PI Loop Filter
    % Ts = 30ms, damping ration = 0.7
    % We get natural frequency = 110, Kp=166.6 and Ki=27755.55
    % B0=166.877556 & B1=-166.322444
    %
    % PI(s) = (B0 + B1*(z^-1))/(1-(z^-1))
    % B0 = (2*Kp + Ki*Ts)/2 | B1 = -(2*Kp - Ki*Ts)/2 
    %---------------------------------------------------------------------%
    
    % PI Coefficients Calculation
    
    Kp    =  166.6;
    Ki    =  27755.55;
    
    PI_B0 =  (2*Kp + Ki*Ts)/2;
    PI_B1 = -(2*Kp - Ki*Ts)/2;
    
    ylf(1)    = 1.0*ylf(2)+ PI_B0*ynotch(1)+ PI_B1*ynotch(2);
    
    % Update Ynotch for future use
    
    ynotch(3) = ynotch(2);
    ynotch(2) = ynotch(1);
    ynotch_buff(n+1)= ynotch(1);
    
    ylf(1)      =  min([ylf(1) 200]);
    ylf(2)      =  ylf(1);
    wo          =  wn + ylf(1);
    
    werror(n+1) = (wo-wn)*0.00318309886;
    
    % Integration process
    Mysin(1) = Mysin(2)+wo*Ts*(Mycos(2));
    Mycos(1) = Mycos(2)-wo*Ts*(Mysin(2));
   
    % Limit the oscillator integrators 
    Mysin(1) = max([Mysin(1) -1.0]);
    Mysin(1) = min([Mysin(1)  1.0]);
    Mycos(1) = max([Mycos(1) -1.0]);
    Mycos(1) = min([Mycos(1)  1.0]);
    Mysin(2) = Mysin(1);
    Mycos(2) = Mycos(1);
    
    % Update the output phase
    theta(1) = theta(2) + wo*Ts;
    
    % Output phase reset condition
    if(Mysin(1)>0 && Mysin(2) <=0)
        theta(1)= -pi;
    end
    
    SinGen(n+1)   = Mycos(1);
    Plot_Var(n+1) = Mysin(1);

    
end

% CASE 1 : Phase Jump at the Mid Point
error = Plot_Var-u;
% CASE 2 : Harmonics
% error=Plot_Var-u1;
% CASE 3: Frequency Variations
% error=Plot_Var-u;
% CASE 4: Amplitude Variations
% error=Plot_Var-u1;

figure;
subplot(3,1,1),plot(t,Plot_Var,'r',t,u,'b'),title('SPLL(red) & Ideal Grid(blue)');
subplot(3,1,2),plot(t,error,'r'),title('Error');
subplot(3,1,3),plot(t,u1,'r',t,Plot_Var,'b'),title('SPLL Out(Blue) & Ideal Grid(Red)');
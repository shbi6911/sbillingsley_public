%Bobby Hodgkinson
%Parameters for the simulink simulation

%% define constants for both rigid and flexible arms
   
%use assignin to write these variables to the workspace
   Kg = 33.3;%9.68*5; %48.4;  % total gear ratio (9.68:1 motor gear box * 5 motor gear to hub gear)
   Km = 0.0401;%0.01065;   %.00767(old);    % 0.01065 default motor constant, V/(rad/s) or Nm/amp.  From new motor data sheet (see Quansermotor.PNG)
   Rm = 19.2;%3.2934;    %2.6(old);   % 3.2934 default output resistance, ohms. Assume this is the 'dynamic resistance' (see Quansermotor.PNG)

%%%%%%%%%%INERTIAS%%%%%%%%%%%
%%%%Assume the same as the ITLL quanser modules.  More testing might be necessary
% Taken from Quanser manual, included base intertia, gears and
% motor inertia.   
J_hub = 0.0005; % base inertia, kgm^2
        
%% define constants for rigid arm
J_load = .0015; %from solidworks load inertia, kgm^2     
J_extra = 0.2*.2794^2; %extra mass
J = J_hub + J_load + J_extra; % effective inertia, kgm^2
B = .05;            %Friction was 0.2


%% define extra constants for flexible arm
L = .45;                    % link length (m), same as the ITLL quanser modules.  
Marm = 0.06;                % link mass of stainless steel ruler
J_arm = (Marm*L^2)/3;       % inertia of the ruler, 0.0041 kgm^2  
Mtip = .050;                % kg. This is the 50 gm tip mass included for stability
J_M = Mtip*L^2;             % 0.01 kgm^2. This is the inertia of the added tip mass
J_L = J_arm + J_M;          % Total load inertia of the flexible link
fc = 1.8;                   % natural frequency, Hz
K_arm = (2*pi*fc)^2*(J_L);  % stiffness (Jl+Jm)
%J_hub is the same as for the rigid arm, they use the same base.

% define constants for simplicity (see equation 22 of ASEN2003_L6_Control.docx)
p1 = -Kg^2*Km^2/(J_hub*Rm); 
q1 = K_arm/(L*J_hub);
r1 = Kg*Km/(J_hub*Rm);
p2 = Kg^2*Km^2*L/(J_hub*Rm);
q2 = -K_arm*(J_hub + J_L)/(J_L*J_hub);
r2 = -Kg*Km*L/(J_hub*Rm);



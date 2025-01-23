function PlotAircraftSim(time, aircraft_state_array, control_input_array, fig, col)

figure(fig(1));
 
subplot(311);
plot(time, aircraft_state_array(:,1),col)
hold on; grid on;
title('Position vs Time')  
ylabel('X [m]')    

subplot(312);
plot(time, aircraft_state_array(:,2),col)
hold on; grid on;
ylabel('Y [m]')    
 
subplot(313);
plot(time, aircraft_state_array(:,3),col)
hold on; grid on;
ylabel('Z [m]')    
xlabel('time [sec]')


figure(fig(2))
 
subplot(311)
plot(time, (180/pi)*aircraft_state_array(:,4),col)
hold on; grid on;
title('Euler Angles vs Time')
ylabel('Roll [deg]')    

subplot(312);
plot(time, (180/pi)*aircraft_state_array(:,5),col)
hold on; grid on;
ylabel('Pitch [deg]')    
 
subplot(313);
plot(time, (180/pi)*aircraft_state_array(:,6),col)
hold on; grid on;
ylabel('Yaw [deg]')    
xlabel('time [sec]')


figure(fig(3));
 
subplot(311);
plot(time, aircraft_state_array(:,7),col)
hold on; grid on;
title('Velocity vs Time')   
ylabel('uE [m/s]')    

subplot(312);
plot(time, aircraft_state_array(:,8),col)
hold on; grid on;
ylabel('vE [m/s]')    
 
subplot(313);
plot(time, aircraft_state_array(:,9),col)
hold on; grid on;
ylabel('wE [m/s]')    
xlabel('time [sec]')


figure(fig(4));
 
subplot(311);
plot(time, (180/pi)*aircraft_state_array(:,10),col)
hold on; grid on;
title('Angular Velocity vs Time');   
ylabel('p [deg/s]')    

subplot(312);
plot(time, (180/pi)*aircraft_state_array(:,11),col)
hold on; grid on;
ylabel('q [deg/s]')    
 
subplot(313);
plot(time, (180/pi)*aircraft_state_array(:,12),col)
hold on; grid on;
ylabel('r [deg/s]')    
xlabel('time [sec]');


figure(fig(5));
hold on; grid on;
plot3(aircraft_state_array(:,1),aircraft_state_array(:,2),-aircraft_state_array(:,3),col)
title('3D Flight Path over Time')
xlabel('X-position (m)');
ylabel('Y-position (m)');
zlabel('Z-position (m)');
%zlim([min(-aircraft_state_array(:,3))-0.005*min(-aircraft_state_array(:,3))...
    %max(-aircraft_state_array(:,3))+ 0.005*max(-aircraft_state_array(:,3))]);

figure(fig(6));

subplot(411);
plot(time, control_input_array(1,:),col)
hold on; grid on;
title('Control Surfaces vs Time') 
ylabel('\delta_e [deg]')    

subplot(412);
plot(time, control_input_array(2,:),col)
hold on; grid on;
ylabel('\delta_a [deg]')      

subplot(413);
plot(time, control_input_array(3,:),col)
hold on; grid on;
ylabel('\delta_r [deg]')       

subplot(414);
plot(time, control_input_array(4,:),col)
hold on; grid on;
ylim([0 1])
ylabel('\delta_t [frac]')     
xlabel('time [sec]');

end





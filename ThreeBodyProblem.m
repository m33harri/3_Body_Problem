%% 3-Body Problem 2D simulation - Tolga Selcuk

% Time For Sim
dt = 0.00005; 

% Mass values for Planets
m1 = 5; m2 = 5; m3 = 9; 
M = [m1, m2, m3]; 

% Initial velocity and relative position conditions
P = [0 25 50; 0 0 0];  % Positions of the bodies
Vo =[-1.47 -1.47 -2; 3.32 3.32 -10.64];     % Initial velocities

% Gravitational constant (normalized to 1 in this simulation)
G = 30;

% Store initial positions for markers
x10 = P(1, 1); y10 = P(2, 1);
x20 = P(1, 2); y20 = P(2, 2);
x30 = P(1, 3); y30 = P(2, 3);

% Set up figure for plotting
figure;
axis([-100 60 -100 60]); % Adjust axis limits as needed
hold on;
grid on; % Add a grid for better visualization
title('Simulation Graphic');
xlabel('x');
ylabel('y');

time = 0:dt:100;

% Preallocate arrays for trajectory plotting
x1 = 0*time; y1 = x1;
x2 = x1;     y2 = x1;
x3 = x1;     y3 = x1;

% Initialize plot handles for bodies
plotHandle1 = plot(P(1,1), P(2,1), 'ro', 'MarkerFaceColor', 'r'); % Body 1
plotHandle2 = plot(P(1,2), P(2,2), 'go', 'MarkerFaceColor', 'g'); % Body 2
plotHandle3 = plot(P(1,3), P(2,3), 'bo', 'MarkerFaceColor', 'b'); % Body 3

tic
for i = 1:1:length(time)

    % Calculate forces on the bodies using pairwise distances
    % Force from body 2 on body 1
    r12 = P(:,2) - P(:,1); % Displacement vector from body 1 to body 2
    distance12 = norm(r12); % Euclidean distance between body 1 and body 2
    F12 = G * M(1) * M(2) * r12 / distance12^3; % Gravitational force from body 2 on body 1

    % Force from body 1 on body 2 (opposite direction)
    F21 = -F12;

    % Force from body 3 on body 1
    r13 = P(:,3) - P(:,1); % Displacement vector from body 1 to body 3
    distance13 = norm(r13); % Euclidean distance between body 1 and body 3
    F13 = G * M(1) * M(3) * r13 / distance13^3; % Gravitational force from body 3 on body 1

    % Force from body 1 on body 3 (opposite direction)
    F31 = -F13;

    % Force from body 3 on body 2
    r23 = P(:,3) - P(:,2); % Displacement vector from body 2 to body 3
    distance23 = norm(r23); % Euclidean distance between body 2 and body 3
    F23 = G * M(2) * M(3) * r23 / distance23^3; % Gravitational force from body 3 on body 2

    % Force from body 2 on body 3 (opposite direction)
    F32 = -F23;

    % Recalculate positions using the net forces
    % Recalculating Position for each body
    P(:,1) = P(:,1) + Vo(:,1)/2 * dt + (F12 + F13) * (dt^2); 
    P(:,2) = P(:,2) + Vo(:,2)/2 * dt + (F21 + F23) * (dt^2); 
    P(:,3) = P(:,3) + Vo(:,3)/2 * dt + (F32 + F31) * (dt^2);
    
    % Recalculating Velocity for each body
    Vo(:,1) = Vo(:,1) + (F12 + F13) * dt;
    Vo(:,2) = Vo(:,2) + (F21 + F23) * dt;
    Vo(:,3) = Vo(:,3) + (F32 + F31) * dt;

    % Store the current positions for trajectory plotting
    x1(i) = P(1,1);
    y1(i) = P(2,1);
    x2(i) = P(1,2);
    y2(i) = P(2,2);
    x3(i) = P(1,3);
    y3(i) = P(2,3);

%    % Update the positions of the bodies in the plot
%    set(plotHandle1, 'XData', P(1,1), 'YData', P(2,1));
%    set(plotHandle2, 'XData', P(1,2), 'YData', P(2,2));
%    set(plotHandle3, 'XData', P(1,3), 'YData', P(2,3));

    % Plot the new positions of the bodies and their trajectories
%    plot(x1(1:t), y1(1:t), 'r', 'LineWidth', 1);
%    plot(x2(1:t), y2(1:t), 'g', 'LineWidth', 1);
%    plot(x3(1:t), y3(1:t), 'b', 'LineWidth', 1);

%    drawnow 
end

toc

%%plot

dtplot = 0.1;

for i=1:floor(dtplot/dt):length(time)
    % Update the positions of the bodies in the plot
   set(plotHandle1, 'XData', x1(i), 'YData', y1(i));
    set(plotHandle2, 'XData', x2(i), 'YData', y2(i));
    set(plotHandle3, 'XData', x3(i), 'YData', y3(i));

    % Plot the new positions of the bodies and their trajectories
    plot(x1(1:i), y1(1:i), 'r');
    plot(x2(1:i), y2(1:i), 'g');
    plot(x3(1:i), y3(1:i), 'b');
     
    title(['Simulation Graphic at t =',num2str(time(i)) ]);
    pause(eps*10^10)
    if (max([abs(x1(i)),abs(x2(i)),abs(x1(i))])>10^5)
        break
    end
end
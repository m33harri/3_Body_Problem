%% 3-Body Problem 2D simulation - Tolga Selcuk

% Time For Sim
dt = 0.01; 

% Mass values for Planets
m1 = 3; m2 = 5; m3 = 10; 
M = [m1, m2, m3]; 

% Initial velocity and relative position conditions
P = [0 -25 50; 0 25 0];  % Positions of the bodies
Vo = zeros(2, 3);     % Initial velocities

% Gravitational constant (normalized to 1 in this simulation)
G = 25;

% Store initial positions for markers
x10 = P(1, 1); y10 = P(2, 1);
x20 = P(1, 2); y20 = P(2, 2);
x30 = P(1, 3); y30 = P(2, 3);

% Set up figure for plotting
figure;
axis([-100 100 -100 100]); % Adjust axis limits as needed
hold on;
grid on; % Add a grid for better visualization
title('Simulation Graphic');
xlabel('x');
ylabel('y');

% Preallocate arrays for trajectory plotting
x1 = []; y1 = [];
x2 = []; y2 = [];
x3 = []; y3 = [];

% Initialize plot handles for bodies
plotHandle1 = plot(P(1,1), P(2,1), 'ro', 'MarkerFaceColor', 'r'); % Body 1
plotHandle2 = plot(P(1,2), P(2,2), 'go', 'MarkerFaceColor', 'g'); % Body 2
plotHandle3 = plot(P(1,3), P(2,3), 'bo', 'MarkerFaceColor', 'b'); % Body 3

for t = 0:dt:30 

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
    x1 = [x1, P(1,1)];
    y1 = [y1, P(2,1)];
    x2 = [x2, P(1,2)];
    y2 = [y2, P(2,2)];
    x3 = [x3, P(1,3)];
    y3 = [y3, P(2,3)];

    % Update the positions of the bodies in the plot
    set(plotHandle1, 'XData', P(1,1), 'YData', P(2,1));
    set(plotHandle2, 'XData', P(1,2), 'YData', P(2,2));
    set(plotHandle3, 'XData', P(1,3), 'YData', P(2,3));

    % Plot the new positions of the bodies and their trajectories
    plot(x1, y1, 'r', 'LineWidth', 1);
    plot(x2, y2, 'g', 'LineWidth', 1);
    plot(x3, y3, 'b', 'LineWidth', 1);

    % Pause for visualization
    pause(0.01); 
end

function [mu, sigma, sigma_points] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model.
% mu: state vector containing robot pose and poses of landmarks obeserved so far
% Current robot pose = mu(1:3)
% Note that the landmark poses in mu are stacked in the order by which they were observed
% sigma: the covariance matrix of the system.
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

% For computing lambda.
global scale;

% Compute sigma points
sigma_points = compute_sigma_points(mu, sigma);

% Dimensionality
n = length(mu);
% lambda
lambda = scale - n;

% TODO: Transform all sigma points according to the odometry command
% Remember to vectorize your operations and normalize angles
% Tip: the function normalize_angle also works on a vector (row) of angles


for i = 1:2*n+1
    theta = sigma_points(3,i);
    %x = sigma_points(1,i);
    %y = sigma_points(2,i);
    %sigma_points(1,i) = x + u.t*cos(theta+u.r1);
    %sigma_points(2,i) = y + u.t*sin(theta+u.r1);
    %sigma_points(3,i) = theta + u.r1 + u.r2;
    addOdo = [u.t*cos(theta+u.r1);u.t*sin(theta+u.r1);u.r1 + u.r2];
    sigma_points(:,i) = sigma_points(:,i) + addOdo;
    normalize_angle(theta);
end
   

% Computing the weights for recovering the mean
wm = [lambda/scale, repmat(1/(2*scale),1,2*n)];
wc = wm;

% TODO: recover mu.
% Be careful when computing the robot's orientation (sum up the sines and
% cosines and recover the 'average' angle via atan2)
summand = zeros(n,1);
summand_x_ba = 0;
summand_y_ba = 0;
mu = zeros(n,1);

for i = 1:2*n+1
    summand_x_ba = summand_x_ba + wm(1,i)*cos(summand(3,i));
    summand_y_ba = summand_y_ba + wm(1,i)*sin(summand(3,i));
end

mu(3,1) = normalize_angle(atan2(summand_y_ba, summand_x_ba));

for i = 1:2*n+1
    summand = summand + wm(1,i)*sigma_points(:,i); 
end
summand(3) = 0;
mu = mu + summand;


% TODO: Recover sigma. Again, normalize the angular difference
sigma = zeros(n);

for i=1:2*n+1
    tmp = sigma_points(:,i) - mu;
    tmp(3) = normalize_angle(tmp(3));
    sigma = sigma + wc(i)* tmp * tmp';
end

% Motion noise
motionNoise = 0.1;
R3 = [motionNoise, 0, 0; 
     0, motionNoise, 0; 
     0, 0, motionNoise/10];
R = zeros(size(sigma,1));
R(1:3,1:3) = R3;

% TODO: Add motion noise to sigma
sigma = sigma + R;

end

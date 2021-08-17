function [mu, sigma] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model,
% mu: 2N+3 x 1 vector representing the state mean
% sigma: 2N+3 x 2N+3 covariance matrix
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

%preparation
t = u.t;
r1 = u.r1;
r2 = u.r2;
theta = mu(3,1);
twoNp3 = size(mu,1);
Fx = [eye(3) zeros(3,twoNp3-3)];
% TODO: Compute the new mu based on the noise-free (odometry-based) motion model
% Remember to normalize theta after the update (hint: use the function normalize_angle available in tools)

mu = mu + Fx * [t*cos(theta+r1);
           t*sin(theta+r1);
           r1+r2];
normalize_angle(theta);

% TODO: Compute the 3x3 Jacobian Gx of the motion model
Gx_low = [1, 0, -t*sin(theta+r1);
          0, 1, t*cos(theta+r1);
          0, 0, 1];


% TODO: Construct the full Jacobian G
G = [Gx_low, zeros(3,twoNp3-3);
     zeros(twoNp3-3,3), eye(twoNp-3)];
 
% Motion noise
motionNoise = 0.1;
R3 = [motionNoise, 0, 0; 
     0, motionNoise, 0; 
     0, 0, motionNoise/10];
R = zeros(size(sigma,1));
R(1:3,1:3) = R3;

% TODO: Compute the predicted sigma after incorporating the motion
sigma = G*sigma*G' + R;

end

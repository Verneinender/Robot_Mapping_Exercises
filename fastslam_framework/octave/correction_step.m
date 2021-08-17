function particles = correction_step(particles, z)

% Weight the particles according to the current map of the particle
% and the landmark observations z.
% z: struct array containing the landmark observations.
% Each observation z(j) has an id z(j).id, a range z(j).range, and a bearing z(j).bearing
% The vector observedLandmarks indicates which landmarks have been observed
% at some point by the robot.

% Number of particles
numParticles = length(particles);

% Number of measurements in this time step
m = size(z, 2);

% TODO: Construct the sensor noise matrix Q_t (2 x 2)
Q_t = 0.1*eye(2);

% process each particle
for i = 1:numParticles
  robot = particles(i).pose;
  % process each measurement
  for j = 1:m
    % Get the id of the landmark corresponding to the j-th observation
    % particles(i).landmarks(l) is the EKF for this landmark
    l = z(j).id;

    % The (2x2) EKF of the landmark is given by
    % its mean particles(i).landmarks(l).mu
    % and by its covariance particles(i).landmarks(l).sigma

    % If the landmark is observed for the first time:
    if (particles(i).landmarks(l).observed == false)

      % TODO: Initialize its position based on the measurement and the current robot pose:
      px = particles(i).pose(1);
      py = particles(i).pose(2);
      ptheta = particles(i).pose(3);
      %lx = px + meas_range * np.cos(ptheta + meas_bearing)
      %ly = py + meas_range * np.sin(ptheta + meas_bearing)
      landmark_x = px + z(j).range * cos(ptheta + z(j).bearing);
      landmark_y = py + z(j).range * sin(ptheta + z(j).bearing);
      particles(i).landmarks(l).mu = [landmark_x; landmark_y];

      % get the Jacobian with respect to the landmark position
      [h, H] = measurement_model(particles(i), z(j));

      % TODO: initialize the EKF for this landmark
      H_inv = pinv(H);
      particles(i).landmarks(l). sigma = H_inv * Q_t * transpose(H_inv);

      % Indicate that this landmark has been observed
      particles(i).landmarks(l).observed = true;

    else

      % get the expected measurement
      [expectedZ, H] = measurement_model(particles(i), z(j));

      % TODO: compute the measurement covariance
      Q = H * particles(i).landmarks(l).sigma * H' + Q_t;

      % TODO: calculate the Kalman gain
      K = particles(i).landmarks(l).sigma * H' * pinv(Q);
      

      % TODO: compute the error between the z and expectedZ (remember to normalize the angle)
      z_diff = [z(j).range; z(j).bearing] - expectedZ;
	    z_diff(2) = normalize_angle(z_diff(2));
      % TODO: update the mean and covariance of the EKF for this landmark
      particles(i).landmarks(l).mu = particles(i).landmarks(l).mu + K * z_diff;
      particles(i).landmarks(l).sigma = particles(i).landmarks(l).sigma - K * H * particles(i).landmarks(l).mu;
      %particles(i).landmarks(l).sigma = (eye(2) - K*H)*particles(i).landmarks(l).sigma;

      % TODO: compute the likelihood of this observation, multiply with the former weight
      % to account for observing several features in one time step
      %fact = 1 / np.sqrt(math.pow(2*math.pi,2) * np.linalg.det(Q))
      %expo = -0.5 * np.dot(delta.T, np.linalg.inv(Q)).dot(delta)
      %weight = fact * np.exp(expo)
      
      temp = 1 / sqrt(det(2 * pi * Q));
      expo = -0.5 * z_diff' * pinv(Q) * z_diff;
      weight = temp * exp(expo);
      particles(i).weight = particles(i).weight * weight;
      

    end

  end % measurement loop
end % particle loop

end

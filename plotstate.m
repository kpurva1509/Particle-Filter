function plotstate(z)
    global determinants;
    global State;
    global Param;
    Sigma = State.Ekf.Sigma;
    mu = State.Ekf.mu;

    robot_pos = mu(1:3);
    robot_variance = Sigma(1:2, 1:2);

    plotcov2d(mu(1), mu(2), Sigma(1:2, 1:2), 'b', 0, 0, 0, 3);

    for index = 1:State.Ekf.nL
        indices = State.Ekf.iL{index};
        map_x = mu(indices(1));
        map_y = mu(indices(2));
        map_sigma = Sigma(indices(1):indices(2), indices(1):indices(2));
        plotcov2d(map_x, map_y, map_sigma, 'r', 0, 0, 0, 3);
        determinants(index,State.Ekf.t) = det(map_sigma);
    end


    title(strcat("Running for DA = ", Param.dataAssociation));
    
    est_x = State.Ekf.mu(1);
    est_y = State.Ekf.mu(2);
    est_theta = State.Ekf.mu(3);
end

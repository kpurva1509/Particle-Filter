function ekfupdate(z)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;

% Check if the identity already exists in the state vector signatures, else addin new landmark
% returns state vector indices pairing observations with landmark
switch lower(Param.dataAssociation)
    case 'known'
        for index = 1:size(z(3, :),2)
            identity = z(3, index);
            exists = false;
            for signature = State.Ekf.sL
                if signature==identity
                    exists = true;
                    break;
                end
            end
            if exists==false
                % Create landmarks and create new indices and return them
                initialize_new_landmark(z(:, index), Param.R);
            end
            
            Li = da_known(z(3,index));
            Q_t = Param.R;
            indices = Li{1};
            robot_x_bar = State.Ekf.mu_bar(1);
            robot_y_bar = State.Ekf.mu_bar(2);
            robot_theta_bar = State.Ekf.mu_bar(3);
            map_x = State.Ekf.mu_bar(indices(1));
            map_y = State.Ekf.mu_bar(indices(2));

            delta_x = map_x - robot_x_bar;
            delta_y = map_y - robot_y_bar;
            delta = [delta_x; delta_y];
            q = delta'*delta;
            z_t_bar = [sqrt(q); minimizedAngle(atan2(delta_y, delta_x) - robot_theta_bar)];
            n_rows = 5;
            n_cols = 3 + State.Ekf.nL*2;
            F_x = zeros(n_rows, n_cols);
            F_x(1:3, 1:3) = eye(3,3);
            F_x(4, indices(1)) = 1;
            F_x(5, indices(2)) = 1;

            H_t = 1/q*[-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0, sqrt(q)*delta_x, sqrt(q)*delta_y;
                delta_y, -delta_x, -q, -delta_y, delta_x]*F_x;
            K_t = State.Ekf.sigma_bar*H_t'*inv(H_t*State.Ekf.sigma_bar*H_t'+Q_t);
            delta_z = z(1:2,index) - z_t_bar;
            delta_z(2) = minimizedAngle(delta_z(2));
            
            State.Ekf.mu_bar = State.Ekf.mu_bar + K_t*(delta_z);
            State.Ekf.mu_bar(3) = minimizedAngle(State.Ekf.mu_bar(3));
            State.Ekf.sigma_bar = (eye(3 + 2*State.Ekf.nL) - K_t*H_t)*State.Ekf.sigma_bar;
        end
        State.Ekf.mu = State.Ekf.mu_bar;
        State.Ekf.Sigma = State.Ekf.sigma_bar;

    case 'known_batch'
        % Update new landmarks before starting batch update
        for index = 1:size(z(3, :), 2)
            identity = z(3, index);
            exists = false;
            for signature = State.Ekf.sL
                if signature==identity
                    exists = true;
                    break;
                end
            end
            if exists==false
                % Create new landmark
                initialize_new_landmark(z(:, index), Param.R);
            end
        end

        % Algorithm for batch update using all landmarks

        % Calculate observations based on the signatures and the predicted states
        z_t = z(1:2, :);
        z_t = z_t(:);
        num_obs = length(z(3, :));
        c_t = z(3, :);
        z_t_bar = [];
        H_t = [];
        Q_t = [];
        robot_x_bar = State.Ekf.mu_bar(1);
        robot_y_bar = State.Ekf.mu_bar(2);
        robot_theta_bar = State.Ekf.mu_bar(3);
        for index = 1:num_obs
            identity = c_t(index);
            Li = da_known(identity);
            indices = Li{1};
            map_x_index = Li{1}(1);
            map_y_index = Li{1}(2);
            map_x = State.Ekf.mu_bar(map_x_index);
            map_y = State.Ekf.mu_bar(map_y_index);
            delta_x = map_x - robot_x_bar;
            delta_y = map_y - robot_y_bar;
            delta = [delta_x; delta_y];
            q = delta'*delta;
            z_t_bar_index = [sqrt(q); minimizedAngle(atan2(delta_y, delta_x) - robot_theta_bar)];
            z_t_bar = [z_t_bar; z_t_bar_index];

            n_rows = 5;
            n_cols = 3 + State.Ekf.nL*2;

            F_x = zeros(n_rows, n_cols);
            F_x(1:3, 1:3) = eye(3,3);
            F_x(4, indices(1)) = 1;
            F_x(5, indices(2)) = 1;
            
            H_t_temp = 1/q*[-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0, sqrt(q)*delta_x, sqrt(q)*delta_y;
                delta_y, -delta_x, -q, -delta_y, delta_x]*F_x;
            H_t = [H_t; H_t_temp];
            
            Q_t = blkdiag(Q_t, Param.R);
        end
        
        S_t = H_t*State.Ekf.sigma_bar*H_t' + Q_t;
        K_t = State.Ekf.sigma_bar*H_t'*inv(S_t);
        State.Ekf.mu = State.Ekf.mu_bar + K_t*(z_t - z_t_bar);
        State.Ekf.Sigma = (eye(3 + 2*State.Ekf.nL) - K_t*H_t)*State.Ekf.sigma_bar;
    case 'nn'
        % Ignore the landmark indices here since we shouldn't use it
        actual_landmarks = z(3, :);
        z  = z(1:2,:);

        % I dont use da_nn function, since it would complicate my logic to implement nn 
        %Li = da_nn(z(1:2,:), Param.R);

        num_obs = size(z, 2);

        landmark_associates = [];

        robot_x_bar = State.Ekf.mu_bar(1);
        robot_y_bar = State.Ekf.mu_bar(2);
        robot_theta_bar = State.Ekf.mu_bar(3);

        for index = 1:num_obs
            z_t = z(1:2, index);

            curr_num_landmarks = State.Ekf.nL;
            % Initialize mahalanobis_distances for each landmarks as well as newly assumed landmark 
            mahalanobis_distances = [];

            for landmark_index = 1:curr_num_landmarks
                indices = State.Ekf.iL{landmark_index};
                map_x_index = indices(1);
                map_y_index = indices(2);

                map_x = State.Ekf.mu_bar(map_x_index);
                map_y = State.Ekf.mu_bar(map_y_index);
                
                % Calculate mahalanobis_distances for each landmarks among my context
                delta_x = map_x - robot_x_bar;
                delta_y = map_y - robot_y_bar;
                delta = [delta_x; delta_y];
                q = delta'*delta;
                z_t_bar = [sqrt(q); minimizedAngle(atan2(delta_y, delta_x) - robot_theta_bar)];

                n_rows = 5;
                n_cols = 3 + State.Ekf.nL*2;
                F_x = zeros(n_rows, n_cols);
                F_x(1:3, 1:3) = eye(3,3);
                F_x(4, indices(1)) = 1;
                F_x(5, indices(2)) = 1;


                H_t = 1/q*[-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0, sqrt(q)*delta_x, sqrt(q)*delta_y;
                    delta_y, -delta_x, -q, -delta_y, delta_x]*F_x;
                C_t = H_t*State.Ekf.sigma_bar*H_t' + Param.R;
                z_delta = z_t - z_t_bar;
                z_delta(2) = minimizedAngle(z_delta(2));
                m_distance = (z_delta)'*inv(C_t)*(z_delta);
                mahalanobis_distances = [mahalanobis_distances, m_distance];
            end
            
            if strcmp(Param.choice, 'sim')
                chi_val = chi2inv(0.985, 9);
            else
                chi_val = 200;
            end
            
            mahalanobis_distances = [mahalanobis_distances, chi_val];

            % Calculate the index of associated landmark among the distances by checking which index has min
            [val, landmark_index] = min(mahalanobis_distances);

            if landmark_index > State.Ekf.nL
                % Add new landmark and associate observation with that one
                initialize_new_landmark([z_t; landmark_index], Param.R)
            end

            robot_x_bar = State.Ekf.mu_bar(1);
            robot_y_bar = State.Ekf.mu_bar(2);
            robot_theta_bar = State.Ekf.mu_bar(3);
            
            State.Ekf.actual_associations(end+1) = actual_landmarks(index);
            try
                State.Ekf.nn_associations(end+1) = Param.landmark_mapping(landmark_index);
            catch
                State.Ekf.nn_associations(end+1) = 0;
            end
            
            Li = State.Ekf.iL{landmark_index};
            
            
            indices = Li;
            Q_t = Param.R;
            map_x = State.Ekf.mu_bar(indices(1));
            map_y = State.Ekf.mu_bar(indices(2));
            delta_x = map_x - robot_x_bar;
            delta_y = map_y - robot_y_bar;
            delta = [delta_x; delta_y];
            q = delta'*delta;
            z_t_bar = [sqrt(q); minimizedAngle(atan2(delta_y, delta_x) - robot_theta_bar)];
            n_rows = 5;
            n_cols = 3 + State.Ekf.nL*2;
            F_x = zeros(n_rows, n_cols);
            F_x(1:3, 1:3) = eye(3,3);
            F_x(4, indices(1)) = 1;
            F_x(5, indices(2)) = 1;

            H_t = 1/q*[-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0, sqrt(q)*delta_x, sqrt(q)*delta_y; delta_y, -delta_x, -q, -delta_y, delta_x]*F_x;
            K_t = State.Ekf.sigma_bar*H_t' / (H_t*State.Ekf.sigma_bar*H_t'+Q_t);
            delta_z = z(1:2,index) - z_t_bar;
            delta_z(2) = minimizedAngle(delta_z(2));
            
            State.Ekf.mu_bar = State.Ekf.mu_bar + K_t*(delta_z);
            State.Ekf.mu_bar(3) = minimizedAngle(State.Ekf.mu_bar(3));
            State.Ekf.sigma_bar = (eye(3 + 2*State.Ekf.nL) - K_t*H_t)*State.Ekf.sigma_bar;
        end
        State.Ekf.mu = State.Ekf.mu_bar;
        State.Ekf.Sigma = State.Ekf.sigma_bar;
        
    case 'jcbb'
        Li = da_jcbb(z(1:2,:), Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end




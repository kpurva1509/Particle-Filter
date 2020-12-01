function Li = da_known(z)
% EKF-SLAM data association with known correspondences

global Param;
global State;

% Since all landmarks are already updated, lets check to associate them appropriately

Li = {};

for identity = z
    signature_index = 1;
    for signature = State.Ekf.sL
        if signature==identity
            Li{end+1} = State.Ekf.iL{signature_index};
        end
        signature_index = signature_index + 1;
    end
end


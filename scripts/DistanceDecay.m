function Q_new = DistanceDecay(Q, Coordinates, lambda, epsilon,Pref)
    % Q: Original transition probability matrix (nearest neighbors)
    % lat, lon: Coordinates of each patch
    % lambda: Decay rate for distance-based movement probability

    M = size(Q, 1); % Number of patches
    Q_new = zeros(M, M); % New movement probability matrix
    lat = Coordinates(:,2);
    lon = Coordinates(:,1);

    % Compute distance matrix D
    D = zeros(M, M);
    for i = 1:M
        for j = 1:M
            D(i, j) = haversine_distance(lat(i), lon(i), lat(j), lon(j));
        end
    end

    % Apply distance scaling while keeping directional probabilities
    for i = 1:M
        Q_new(i,:)=Q(i,:);
        % Compute new probabilities with exponential decay
        long_dist_probs = exp(-lambda * D(i, :));  
        long_dist_probs(Q(i,:)>0)=0;
        long_dist_probs(i) = 0;
        if all(long_dist_probs==0)
            long_dist_probs = 0;
        else
        long_dist_probs = epsilon*(long_dist_probs/sum(long_dist_probs));
        end
        
        % Normalize probabilities to sum to 1
        Q_new(i,:) = Q_new(i,:) + long_dist_probs;
        Q_new(i,:) = Q_new(i,:)/sum(Q_new(i,:));
    end
end

% Haversine distance function for real-world distances (km)
function d = haversine_distance(lat1, lon1, lat2, lon2)
    R = 6371; % Earth radius in km
    dlat = deg2rad(lat2 - lat1);
    dlon = deg2rad(lon2 - lon1);
    a = sin(dlat/2).^2 + cos(deg2rad(lat1)) .* cos(deg2rad(lat2)) .* sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    d = R * c;
end
function [A_eps, A_BN_tilde] = star_sensor(stars, A_BN, FOV)

    scale = 1; % scale of the celestial sphere 
    %n_sensor_B = n_sensor_B./norm(n_sensor_B);
    FOV_rad = deg2rad(FOV);
    toll_FOV = FOV_rad*4/100;

    %% sensor geometry
    n_sens_B = [0; 0; 1];
    n_sens_N = A_BN' * n_sens_B; 
    x_sens_N = A_BN' * [1; 0; 0]; % for sake of simplicity directed as x_body
    y_sens_N = A_BN' * [0; 1; 0];

    max_mag = max(stars.mag);
    min_mag = min(stars.mag);
    f = .5*scale ./ tan(FOV_rad/2); % [m]

    %% print the celestial sphere with the stars
    N_stars = length(stars.mag(stars.mag >= min_mag));
    fprintf("Using %d stars from catalogue. [magitude from %.3f to %.3f]\n", N_stars, min_mag, max_mag);

    filtered_catalogue.ra  = stars.ra(find(stars.mag >= min_mag));
    filtered_catalogue.dec = stars.dec(find(stars.mag >= min_mag));
    filtered_catalogue.mag = stars.mag(stars.mag >= min_mag);


    %% plot the celestial sphere IN INERTIAL EARTH CENTRE
    theta_vect = linspace(0,2*pi);
    phi_vect = linspace(-pi,pi);

    [theta_mesh, phi_mesh] = meshgrid(theta_vect,phi_vect);
    [xHSph,yHSph,zHSph] = sph2cart(theta_mesh,phi_mesh,scale);

    figure(1)
    subplot(1,2,1)
    hSphAx = surf(xHSph,yHSph,zHSph,'FaceColor',[0 1/3 1],'EdgeColor','none');
    alpha(hSphAx,.7);
    hold on;
    axis image;
    title("Celestial Sphere with "+num2str(N_stars)+" stars");

    plot3(0,0,0,'k+','LineWidth',1);

    right_ascension = deg2rad(360 - filtered_catalogue.ra);
    declination     = deg2rad(filtered_catalogue.dec);
    [filtered_catalogue.x,filtered_catalogue.y,filtered_catalogue.z] = sph2cart(right_ascension,declination,scale);
    plot3(filtered_catalogue.x,filtered_catalogue.y,filtered_catalogue.z,'.g','LineWidth',0.1);


    %% transform the pointing vector of sensor into RA & DEC
    n_sens_N_ra  = rad2deg(atan2(n_sens_N(2), n_sens_N(1)));
    n_sens_N_dec = rad2deg(atan2(n_sens_N(3), sqrt(n_sens_N(2)^2+n_sens_N(1)^2)));
    


    %% plot the FOV cone
    [X,Y,Z] = cone(FOV, n_sens_N, scale);
    surf(X,Y,Z,'FaceColor','#63666A','LineStyle','none');
    plot3(n_sens_N(1)*scale,n_sens_N(2)*scale,n_sens_N(3)*scale,'+r','LineWidth',2);

    xlabel("X");
    ylabel("Y");
    zlabel("Z");

    %% plot the image on the sensor
    vectorial_stars = [filtered_catalogue.x(:),filtered_catalogue.y(:),filtered_catalogue.z(:)];
    vectorial_stars = vectorial_stars./vecnorm(vectorial_stars')';
    visible_stars.id = find( acos( (vectorial_stars*n_sens_N(:)) ) <= FOV_rad-toll_FOV ); % sistema con toll

    fprintf("The sensor has %d stars in its FOV.\n", max(size(visible_stars.id)));

    visible_stars.ra  = filtered_catalogue.ra(visible_stars.id);
    visible_stars.dec = filtered_catalogue.dec(visible_stars.id);
    visible_stars.x   = filtered_catalogue.x(visible_stars.id);
    visible_stars.y   = filtered_catalogue.y(visible_stars.id);
    visible_stars.z   = filtered_catalogue.z(visible_stars.id);

    plot3(visible_stars.x,visible_stars.y,visible_stars.z,'or','LineWidth',0.5);

    for i = 1:length(visible_stars.id)
        obs_star_B.x(i) = A_BN(1,:) * [visible_stars.x(i),visible_stars.y(i),visible_stars.z(i)]';
        obs_star_B.y(i) = A_BN(2,:) * [visible_stars.x(i),visible_stars.y(i),visible_stars.z(i)]';
        obs_star_B.z(i) = A_BN(3,:) * [visible_stars.x(i),visible_stars.y(i),visible_stars.z(i)]';
    
        %plot3(sensor_star.x,sensor_star.y,sensor_star.z,'^c','LineWidth',0.5);
    end
    
    lambda = asin(-obs_star_B.y./scale);
    phi    = atan2(-obs_star_B.x,obs_star_B.z);

    u = f.*tan(phi);
    v = f.*tan(lambda)./cos(phi);
    
    theta = linspace(0,2*pi);
    
    subplot(1,2,2)
    plot(scale.*cos(theta),scale.*sin(theta), 'r-','LineWidth',0.4);
    hold on;
    plot(0,0 ,'r+','LineWidth',1);
    xline(0);
    yline(0);
    axis image;
    plot(u, v, '*g');
    title("Stars on Sensor: "+num2str(length(visible_stars.id))+", Ideal Stars VS Real Stars");
    xlim([-scale, scale]);
    ylim([-scale, scale]);
    xlabel("u [-]");
    ylabel("v [-]")

    %% ADDING ERROR
    % lens error (systematic error) - remove if interested in A_eps
        const_u = 1e-3 * scale;
        const_v = 1e-4 * scale;
    % noise (random error)
        sigma = .007;

    % total error
    u_real = normrnd(u, sigma) + const_u;
    v_real = normrnd(v, sigma) + const_v;


    figure(1)
    subplot(1,2,2)
        plot(u_real, v_real, '.b');

    %% ATTITUDE DETERMINATION
    phi_s = atan(u_real./f);
    lambda_s = atan(v_real./f .* cos(phi_s));

    O_hat = [-sin(phi_s(:)).*cos(lambda_s(:)), -sin(lambda_s(:)), cos(phi_s(:)).*cos(lambda_s(:))]';

    sensor_stars = [visible_stars.x(:)'; visible_stars.y(:)'; visible_stars.z(:)'];
    
    A_BN_tilde = O_hat * ( pinv(sensor_stars) ); %% using pseudoinverse
    
    % orthnormalise A_BN_tilde
        
    A_BN
    A_BN_tilde
    
    n_sens_N
    n_sens_N_tilde = A_BN_tilde' * n_sens_B
    
    A_eps = A_BN_tilde*(A_BN')

    figure(1)
    subplot(1,2,1)
    plot3(n_sens_N_tilde(1),n_sens_N_tilde(2),n_sens_N_tilde(3),'xk','LineWidth',2);



%     %% USING STATISTICAL METHODS:
% 
%     %alpha_param = 1;
%     B = O_hat*sensor_stars';
%     [U,~] = eig(B*B');
%     [V,~] = eig(B'*B);
%     D = diag([1 1 det(U)*det(V)]);
% 
%     A_BN_tilde_stat = U*D*V'
%     n_sens_N_tilde_stat = A_BN_tilde_stat' * n_sens_B
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X3,Y3,Z3] = cone(theta,dir,h)
%function: cone. Generates a right circular cone with axis orientation
%aperture angle and height specified.
% inputs : theta: aperture angle of cone
%          d    : vector describing orientation of axis of cone.
%          h    : height of cone 
    
    r = h*tan(pi*theta/180);
    m = h/r;
    [R,A] = meshgrid(linspace(0,r,11),linspace(0,2*pi,41));
    % Generate cone about Z axis with given aperture angle and height
    X = R .* cos(A);
    Y = R .* sin(A);
    Z = m*R;
    % Cone around the z-axis, point at the origin
    % find coefficients of the axis vector xi + yj + zk
    x = dir(1);
    y = dir(2);
    z = dir(3);
    
    % find angle made by axis vector with X axis
    phix = atan2(y,x);
    % find angle made by axis vector with Z axis
    phiz = atan2(sqrt(x^2 + y^2),(z));
    
    % Rotate once about Z axis 
    X1 = X*cos(phiz)+Z*sin(phiz);
    Y1 = Y;
    Z1 = -X*sin(phiz)+Z*cos(phiz);
    % Rotate about X axis
    X3 = X1*cos(phix)-Y1*sin(phix);
    Y3 = X1*sin(phix)+Y1*cos(phix);
    Z3 = Z1;
    
end
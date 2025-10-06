clear; close all;

section = input('Enter the question no. you want the answer for: ');

if ismember(section, [1,2,3,4])
    A = [0 1 0 0; 0 0 1 0; 0 0 0 1; -1 -5 -6 -4];
    B = [0; 0; 0; 1];
    C = [1 0 0 0];
    D = 0;
    sys = ss(A, B, C, D);
end

%% Question-1
if section == 1
    Wc_1 = lyap(A, B*B');
    fprintf('Rank of Controllability matrix via solving Lyapunov Eqn: %d \n', rank(Wc_1))
    Wc_2 = gram(sys, 'c');
    fprintf('Rank of Controllability matrix via "gram" command: %d \n', rank(Wc_2))
    disp('Difference b/w both Controllability matrices')
    disp(abs(Wc_1 - Wc_2))

    Wo_1 = lyap(A', C'*C);
    fprintf('Rank of Observability matrix via solving Lyapunov Eqn: %d \n', rank(Wo_1))
    Wo_2 = gram(sys, 'o');
    fprintf('Rank of Observability matrix via "gram" command: %d \n', rank(Wo_2))
    disp('Difference b/w both Observability matrices')
    disp(abs(Wo_1 - Wo_2))
end
%% Question-2
if section == 2
    Wc = gram(sys, 'c');
    Wo = gram(sys, 'o');

    eig_c = eig(Wc);
    eig_o = eig(Wo);

    figure
    plot(1:length(eig_c), eig_c, '-o', 'DisplayName', 'Eigen values of Wc', 'LineWidth', 1.5)
    hold on
    plot(1:length(eig_o), eig_o, '-o', 'DisplayName', 'Eigen values of Wo', 'LineWidth', 1.5)
    xlabel('Eigen Value No.')
    ylabel('Eigen Values')
    title('Eigen Values of original grammians')
    legend
    grid on

    T = diag([0.1 0.2 0.3 0.5]);
    Wc_trans = inv(T) * Wc * inv(T');
    Wo_trans = T' * Wo * T;

    eig_c_trans = eig(Wc_trans);
    eig_o_trans = eig(Wo_trans);

    figure
    plot(1:length(eig_c_trans), eig_c_trans, '-o', 'DisplayName', 'Eigen values of Wc', 'LineWidth', 1.5)
    hold on
    plot(1:length(eig_o_trans), eig_o_trans, '-o', 'DisplayName', 'Eigen values of Wo', 'LineWidth', 1.5)
    xlabel('Eigen Value No.')
    ylabel('Eigen Values')
    title('Eigen Values of transformed grammians')
    legend
    grid on
end
%% Question-3
if section == 3 || section == 4
    Wc = gram(sys, 'c');
    Wo = gram(sys, 'o');

    Rc = chol(Wc, 'lower');
    Ro = chol(Wo, 'lower');

    [U, S, V] = svd(Ro' * Rc);
    T1 = Rc * V * inv(S)^0.5;
    T1_inv = inv(S)^0.5 * U' * Ro';

    A_bal = T1_inv * A * T1;
    B_bal = T1_inv * B;
    C_bal = C * T1;
    sys_bal = ss(A_bal, B_bal, C_bal, D);

    Wc_1 = gram(sys_bal,'c');
    Wo_1 = gram(sys_bal,'o');
    hsv1 = sqrt(eig(Wc_1 * Wo_1));

    [sysb, hsv2, T2, T2_inv] = balreal(sys);

    if section == 3
        fprintf('T and T_inv  obtained from manual balanced realization are \n')
        disp('T = ');
        disp(T1);
        disp('T_inv = ');
        disp(T1_inv);

        fprintf('T and T_inv  obtained from "balreal" command are \n')
        disp('T = ');
        disp(T2);
        disp('T_inv  = ');
        disp(T2_inv);

        figure
        plot(1:length(hsv1), hsv1, 'r-o', 'DisplayName', 'HSVs of manual balanced realization', 'LineWidth', 1.5)
        hold on
        plot(1:length(hsv2), hsv2, 'g', 'DisplayName', 'HSVs of balreal command', 'LineWidth', 1.5)
        xlabel('Index')
        ylabel('Value')
        %title('Comparison of HSVs from manual balanced realization and balreal command')
        legend
        grid on

    else
        r = 2;

        A_trunc = A_bal(1:r, 1:r);
        B_trunc = B_bal(1:r, :);
        C_trunc = C_bal(:, 1:r);
        D_trunc = D;
        sys_baltrunc = ss(A_trunc, B_trunc, C_trunc, D_trunc);

        sys_btrunc = balred(sysb, r);

        [y1,t1] = step(sys_baltrunc, 20); % manual balanced truncation
        [y2,t2] = step(sys_btrunc, 20);   % balred truncation

        figure;
        step(sys, 'b', 20); % 20 seconds
        hold on
        plot(t1, y1, 'r-o', 'LineWidth', 1.5, 'MarkerIndices', 1:10:length(t1));
        plot(t2, y2, 'g', 'LineWidth', 1.5);
        legend('Full system', 'Manual truncation order 2', ...
            'balred truncation order 2')
        xlabel('Time')
        ylabel('Output')
        title('Step response comparison')
        hold off
        grid on

        lower_bound = hsv2(r+1);
        fprintf('Lower Bound for error: %f\n', lower_bound);
        upper_bound = 2 * sum(hsv2(r+1:end));
        fprintf('Upper Bound for error: %f\n', upper_bound);

        Error1 = norm(sys - sys_baltrunc, inf); % manual balanced truncation
        Error2 = norm(sys - sys_btrunc, inf);   % balred truncation

        fprintf('Actual H-infinity norm of error: %f\n', Error1);
        fprintf('Actual H-infinity norm of error: %f\n', Error2);
    end
end
%% Question-5
if section == 5
    load ./'mat files'/mass_spring_damper_matrices.mat

    eig_A = eig(A);
    if real(eig_A) > 0
        disp('Modes are unstable')
    else
        disp('All the modes are stable')
    end

    sys = ss(A, B, C, D);
    Wc = gram(sys, 'c');
    Wo = gram(sys, 'o');

    order = 12;
    if order == rank(ctrb(A, B))
        disp('System is controllable')
    else
        disp('System is not controllable')
    end

    if order == rank(obsv(A, C))
        disp('System is observable')
    else
        disp('System is not observable')
    end

    [sysb, hsv] = balreal(sys);
    figure
    plot(1:length(hsv), hsv, 'r-o', 'LineWidth', 1.5)
    xlabel('Index')
    ylabel('Value')
    title('Hankel Singular Values of the system')
    grid on

    r = 2; % r = 4 provides the same result due to pole-zero cancellation
    sys_br = balred(sysb, r);
    [y1,t1] = step(sys, 5);
    [y2,t2] = step(sys_br, 5);

    rel_Hinf_err = norm(sys - sys_br, inf)/norm(sys, inf);
    fprintf('Relative H_inf error between the two systems is: %d \n', rel_Hinf_err)

    % Manual Balanced Truncation
    sys_min = minreal(sys);  % Let's remove the uncontrollable states for now
    Wc_1 = lyap(sys_min.A, sys_min.B*(sys_min.B)');
    Wo_1 = lyap((sys_min.A)', (sys_min.C)'*sys_min.C);

    Rc = chol(Wc_1, 'lower');
    Ro = chol(Wo_1, 'lower');

    [U, S, V] = svd(Ro' * Rc);
    T1 = Rc * V * inv(S)^0.5;
    T1_inv = inv(S)^0.5 * U' * Ro';

    A_bal = T1_inv * sys_min.A * T1;
    B_bal = T1_inv * sys_min.B;
    C_bal = sys_min.C * T1;

    r = 2;
    Ar = A_bal(1:r, 1:r);
    Br = B_bal(1:r, :);
    Cr = C_bal(:, 1:r);
    Dr = sys_min.D;
    sys_balr = ss(Ar, Br, Cr, Dr);

    [y3,t3] = step(sys_balr, 5);

    figure;
    for i = 1:6
        subplot(3,2,i)
        plot(t1, y1(:,i), 'r-*', 'LineWidth', 1.2, ...
            'MarkerIndices', 1:10:length(t1));
        hold on
        plot(t2, y2(:,i), 'b-s', 'LineWidth', 1.2, ...
            'MarkerIndices', 1:10:length(t1));
        plot(t3, y3(:,i), 'g', 'LineWidth', 1.5);
        grid on

        xlabel('Time (s)')
        ylabel(['y' num2str(i)])
        legend('full order system','balred truncation','manual truncation','Location','best')
    end
end
%% Question-5
if section == 6
    load ./'mat files'/cylinder_flow_snapshots.mat

    [U, S, V] = svd(X, 'econ');
    sing_vals = diag(S);
    explained_var = (sing_vals.^2) / sum(sing_vals.^2); % Explained variance
    cumulative_var = cumsum(explained_var) * 100;       % Cumulative explained variance
    r = find(cumulative_var > 99.9, 1);
    fprintf('r = %d for capturing 99.9%% variance\n', r);

    Ur = U(:,1:r); Sr = S(1:r,1:r); Vr = V(:,1:r);
    A_hat = Ur' * Y * Vr * inv(Sr);

    [W, L] = eig(A_hat);
    phi = Y * Vr * inv(Sr) * W;

    lambda = diag(L);
    % Assuming dt = 1
    omega = log(lambda); % continuous-time eigenvalues
    b = phi \ X(:,1);    % initial condition coefficients

    time = 0:(size(X,2)-1);
    X_hat = zeros(size(X));
    for j = 1:length(time)
        X_hat(:,j) = phi * (diag(exp(omega*time(j))) * b); % snapshot at time t(k)
    end

    vortmin = -5; vortmax = 5;
    plotCylinder(X, X_hat, nx, ny, r, vortmin, vortmax)

    % Errors w.r.t full order models
    r_list = 10:10:150;
    errors = zeros(size(r_list));
    for k = 1:length(r_list)
        r = r_list(k);

        Ur = U(:,1:r);
        Sr = S(1:r,1:r);
        Vr = V(:,1:r);

        A_hat = Ur' * Y * Vr * inv(Sr);
        [W, L] = eig(A_hat);
        phi = Y * Vr * inv(Sr) * W;

        lambda = diag(L);
        omega = log(lambda);
        b = phi \ X(:,1);

        % Reconstruct snapshots
        time = 0:(size(X,2)-1);
        X_hat = zeros(size(X));
        for j = 1:length(time)
            X_hat(:,j) = phi * (diag(exp(omega*time(j))) * b);
        end

        % Relative frobenius norm error
        errors(k) = norm(X - X_hat, 'fro') / norm(X, 'fro');
    end

    figure;
    semilogy(r_list, errors, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8)
    xlabel('Reduced-order model size (r)')
    ylabel('Relative Frobenius error (log scale)')
    title('Model reduction error vs truncation rank')
    grid on
end
%% Question-7
if section == 7
    % Parameters
    m = 1.0; k = 100; c = 0.5;
    dt = 0.01;
    msnap = 200;
    x0 = [1;0];
    sigma = 0.02;  % measurement noise
    NMC = 500;     % Monte Carlo trials
    r = 2;         % rank for TLS-DMD

    % Continuous and discrete system
    Ac = [0 1; -k/m -c/m];
    Ad = expm(Ac*dt);
    lambda_true = sort(eig(Ad));

    % Generating snapshots
    X = zeros(2, msnap-1);
    Y = zeros(2, msnap-1);
    xk = x0;
    for k_snap = 1:msnap-1
        X(:,k_snap) = xk;
        xk = Ad * xk;
        Y(:,k_snap) = xk;
    end

    eig_DMD_all = zeros(r,NMC);
    eig_TLS_all = zeros(r,NMC);

    % Monte Carlo trials
    for trial = 1:NMC
        % Adding Gaussian noise to snapshots
        E1 = sigma*randn(size(X));
        E2 = sigma*randn(size(Y));
        X_noisy = X + E1;
        Y_noisy = Y + E2;

        % Standard DMD
        [Ux, Sx, Vx] = svd(X_noisy,'econ');
        Ur = Ux(:,1:r); Sr = Sx(1:r,1:r); Vr = Vx(:,1:r);
        Atilde = Ur' * Y_noisy * Vr / Sr; % projected DMD
        eig_DMD_all(:,trial) = sort(eig(Atilde));

        % TLS-DMD
        Z = [X_noisy; Y_noisy];
        [Uz,Sz,Vz] = svd(Z,'econ');
        Uz_r = Uz(:,1:r); Vz_r = Vz(:,1:r); Sz_r = Sz(1:r,1:r);

        % Partition Uz_r: top = U11, bottom = U21
        U11 = Uz_r(1:size(X,1), :);
        U21 = Uz_r(size(X,1)+1:end, :);

        Atls = U21 / U11;   % TLS-DMD operator
        eig_TLS_all(:,trial) = sort(eig(Atls));
    end

    mean_DMD = mean(eig_DMD_all,2);
    std_DMD  = std(eig_DMD_all,0,2);
    mean_TLS = mean(eig_TLS_all,2);
    std_TLS  = std(eig_TLS_all,0,2);

    figure; hold on
    % True eigenvalues
    plot(real(lambda_true), imag(lambda_true), 'k*', 'MarkerSize',10, 'DisplayName','True eigenvalues');
    % Standard DMD
    errorbar(real(mean_DMD), imag(mean_DMD), imag(std_DMD), imag(std_DMD), real(std_DMD), real(std_DMD), ...
        'o','LineWidth',1.2,'DisplayName','Standard DMD');
    % TLS-DMD
    errorbar(real(mean_TLS), imag(mean_TLS), imag(std_TLS), imag(std_TLS), real(std_TLS), real(std_TLS), ...
        's','LineWidth',1.2,'DisplayName','TLS-DMD');

    xlabel('Real part'); ylabel('Imag part')
    title('Eig values of a noisy SDOF system')
    legend
    grid on
    axis equal

    figure('Position', [100, 100, 1200, 600]);
    % 1st eigenvalue
    subplot(1,2,1); hold on
    plot(real(lambda_true(1)), imag(lambda_true(1)), 'k*', 'MarkerSize',10, 'LineWidth',1.5, 'DisplayName','True eig 1');
    scatter(real(eig_DMD_all(1,:)), imag(eig_DMD_all(1,:)), 25, [1 0.4 0.4], 'filled', 'MarkerFaceAlpha',0.15, 'DisplayName','Std DMD');
    scatter(real(eig_TLS_all(1,:)), imag(eig_TLS_all(1,:)), 25, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha',0.15, 'DisplayName','TLS-DMD');
    plot(real(mean_DMD(1)), imag(mean_DMD(1)), 'ro', 'MarkerSize',8, 'LineWidth',2, 'DisplayName','Std DMD mean');
    plot(real(mean_TLS(1)), imag(mean_TLS(1)), 'bs', 'MarkerSize',8, 'LineWidth',2, 'DisplayName','TLS-DMD mean');
    xlabel('Real part'); ylabel('Imag part'); title('Eigenvalue 1')
    legend; grid on; axis equal;

    % 2nd eigenvalue
    subplot(1,2,2); hold on
    plot(real(lambda_true(2)), imag(lambda_true(2)), 'k*', 'MarkerSize',10, 'LineWidth',1.5, 'DisplayName','True eig 2');
    scatter(real(eig_DMD_all(2,:)), imag(eig_DMD_all(2,:)), 25, [1 0.4 0.4], 'filled', 'MarkerFaceAlpha',0.15, 'DisplayName','Std DMD');
    scatter(real(eig_TLS_all(2,:)), imag(eig_TLS_all(2,:)), 25, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha',0.15, 'DisplayName','TLS-DMD');
    plot(real(mean_DMD(2)), imag(mean_DMD(2)), 'ro', 'MarkerSize',8, 'LineWidth',2, 'DisplayName','Std DMD mean');
    plot(real(mean_TLS(2)), imag(mean_TLS(2)), 'bs', 'MarkerSize',8, 'LineWidth',2, 'DisplayName','TLS-DMD mean');
    xlabel('Real part'); ylabel('Imag part'); title('Eigenvalue 2')
    legend; grid on; axis equal;

    sgtitle(sprintf('For measurement noise = %.e', sigma));
end
%% Question-8
if section == 8
    % Parameters
    n = 101;          % Number of grid points
    x = linspace(0,1,n)';
    dx = x(2)-x(1);
    nu = 0.1;         % Viscosity
    dt = 0.0001;      % Time step
    Tfinal = 1;       % Final time
    nt = round(Tfinal/dt);
    u0 = sin(pi*x);   % Initial condition

    e = ones(n,1);    % derivative matrices (central difference)
    D1 = spdiags([-0.5*e, 0*e, 0.5*e], -1:1, n, n)/dx;  % First derivative
    D1(1,:) = 0; D1(end,:) = 0;                         % Dirichlet BC
    D2 = spdiags([e, -2*e, e], -1:1, n, n)/(dx^2);      % Second derivative
    D2(1,:) = 0; D2(end,:) = 0;                         % Dirichlet BC

    % Time integration (Forward Euler)
    U = zeros(n, nt+1);
    U(:,1) = u0;
    for i = 1:nt
        u = U(:,i);
        dudt = - u .* (D1*u) + nu*(D2*u);
        U(:, i+1) = u + dt*dudt;
    end

    % POD Approximation
    snapshots = U;  % Snapshot matrix for POD
    [Phi, S, ~] = svd(snapshots, 'econ');

    singular_vals = diag(S);
    energy = cumsum(singular_vals.^2)/sum(singular_vals.^2);
    r_best = find(energy >= 0.999, 1);  % Select r modes capturing 99.9% energy
    fprintf(['To capture 99.9%% energy of the original system, we need to ' ...
        'atleast truncate it to order %d\n'], r_best)

    r_values = [2, 3, 4, 5, 8, 10];
    colors = lines(length(r_values));
    markers = {'o','s','^','d','v','>','<','p','h','x','*'};
    avg_errors = zeros(size(r_values));

    figure; % Final solutions
    plot(x, U(:,end), 'k-', 'LineWidth', 2); hold on;

    figure; % Error over time
    hold on;
    for idx = 1:length(r_values)
        r = r_values(idx);
        Phi_r = Phi(:,1:r);

        % Reduced diffusion operator
        A = nu * Phi_r' * D2 * Phi_r;

        % Reduced convection operator
        D1Phi = D1 * Phi_r;
        N = zeros(r, r^2);
        for i = 1:r
            for j = 1:r
                col = Phi_r(:,i) .* D1Phi(:,j);
                N(:, (i-1)*r + j) = - Phi_r' * col;
            end
        end

        % ROM Simulation
        u_hat0 = Phi_r' * u0;
        u_hat = zeros(r, nt+1);
        u_hat(:,1) = u_hat0;
        for k = 1:nt
            u_hat_vec = u_hat(:,k);
            u_hat_kron = kron(u_hat_vec, u_hat_vec);
            dudt_hat = A*u_hat_vec + N*u_hat_kron;
            u_hat(:, k+1) = u_hat_vec + dt*dudt_hat;
        end
        U_rom = Phi_r * u_hat;

        % Error analysis
        error = vecnorm(U - U_rom, 2, 1) ./ vecnorm(U,2,1);
        avg_errors(idx) = mean(error);

        % Plots
        figure(1) % Final solution
        plot(x, U_rom(:,end), ...
            'Color', colors(idx,:), ...
            'LineWidth', 1.5, ...
            'Marker', markers{mod(idx-1,length(markers))+1}, ...
            'MarkerIndices', 1:5:n);

        figure(2) % Error over time
        plot(0:dt:Tfinal, error, ...
            'Color', colors(idx,:), ...
            'LineWidth', 1.5, ...
            'Marker', markers{mod(idx-1,length(markers))+1}, ...
            'MarkerIndices', 1:round(nt/20):nt);
    end

    figure(1)
    legend('FOM');
    title('FOM vs ROM Solutions at Final Time');
    legend(arrayfun(@(r) sprintf('ROM (r=%d)',r), r_values, 'UniformOutput',false));
    xlabel('x'); ylabel('u(x,T)');

    figure(2)
    legend(arrayfun(@(r) sprintf('ROM (r=%d)',r), r_values, 'UniformOutput',false));
    xlabel('Time'); ylabel('Relative L2 Error');
    title('Error of ROM vs FOM for different r');

    % Efficiency plot (Accuracy vs Efficiency)
    figure;
    plot(r_values, avg_errors, 'o-', 'LineWidth', 2);
    xlabel('ROM dimension r');
    ylabel('Average Relative L2 Error');
    title('Accuracy vs Efficiency: Error vs ROM size');
    grid on;
end



















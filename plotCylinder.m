% This file is taken from the book Dynamic Mode Decomposition by Nathan Kutz et al. 
% The file was modified by Aniketh Kalur (2020)
% The file was again modified by Charukhesh (2025)

function plotCylinder(X, X_hat, nx, ny, r, vortmin, vortmax)
% plotCylinder - Animates original and reconstructed vorticity fields with a cylinder

% Syntax: plotCylinder(X, X_hat, nx, ny, r, vortmin, vortmax)
% Inputs:
%   X        - Original snapshot matrix (size: nx*ny x time steps)
%   X_hat    - Reconstructed snapshot matrix (same size as X)
%   nx, ny   - Grid dimensions
%   r        - POD/DMD rank used for reconstruction (for title)
%   vortmin  - Minimum vorticity to clip
%   vortmax  - Maximum vorticity to clip

    if nargin < 6
        vortmin = -5; 
        vortmax = 5;
    end

    % Precompute cylinder coordinates
    theta = linspace(0, 2*pi, 200);
    xc = 49 + 25*sin(theta);
    yc = 99 + 25*cos(theta);

    % Initial snapshot
    VORT0 = reshape(real(X(:,1)), nx, ny);
    VORT0(VORT0 > vortmax) = vortmax;
    VORT0(VORT0 < vortmin) = vortmin;

    figure;
    % Original snapshot subplot
    subplot(1,2,1)
    h1 = imagesc(VORT0, [vortmin vortmax]);  
    axis equal; axis off; hold on;
    fill(xc, yc, [.3 .3 .3]);
    plot(xc, yc, 'k', 'LineWidth', 1.2);
    title('Original')

    % Reconstructed snapshot subplot
    subplot(1,2,2)
    h2 = imagesc(VORT0, [vortmin vortmax]);
    axis equal; axis off; hold on;
    fill(xc, yc, [.3 .3 .3]);
    plot(xc, yc, 'k', 'LineWidth', 1.2);
    title(['DMD reconstructed (r = ', num2str(r), ')'])

    colormap(jet)

    % Animation loop
    nFrames = min(size(X,2), 150);  % max 150 frames
    for t = 1:nFrames
        V1 = reshape(real(X(:,t)), nx, ny);
        V2 = reshape(real(X_hat(:,t)), nx, ny);

        % Clip vorticity
        V1(V1 > vortmax) = vortmax; V1(V1 < vortmin) = vortmin;
        V2(V2 > vortmax) = vortmax; V2(V2 < vortmin) = vortmin;

        % Update images
        set(h1, 'CData', V1);
        set(h2, 'CData', V2);

        % Update titles
        subplot(1,2,1), title(['Original, t = ', num2str(t)])
        subplot(1,2,2), title(['DMD reconstructed (r = ', num2str(r), '), t = ', num2str(t)])

        drawnow
        pause(0.01)
    end
end

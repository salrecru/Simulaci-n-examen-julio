% --- Parámetros del sistema ---
    g = 9.81;       % Aceleración de la gravedad (m/s^2)
    a = 0.5;        % Parámetro de la parábola z = a*r^2 (1/m)
    m = 1.0;        % Masa (kg)
    
    % --- CONDICIONES INICIALES ---
    r = 2.0;       % Distancia radial inicial (m)
    vr = 0.0;      % Velocidad radial inicial (m/s)
    phi = 3.0;     % Ángulo inicial (rad)
    phidot = 2.0;    % Velocidad angular inicial (0 para ver la Normal, >0 para giro)
    
    p_phi = m * (r^2) * phidot;
    
    t_final = 15;   % Tiempo total de simulación (s)
    tspan = [0 t_final];
    
    % --- Ecuaciones Diferenciales del Movimiento (Apartado c) ---
    % Y(1) = r, Y(2) = dr/dt, Y(3) = phi
    if phidot == 0
        % Caso Péndulo Plano (phidot = 0): r puede ser positivo y negativo (x con signo)
        ode_fun = @(t, Y) [ ...
            Y(2); ...
            (- 4*a^2*Y(1)*(Y(2)^2) - 2*g*a*Y(1)) / (1 + 4*a^2*Y(1)^2); ...
            0 ...
        ];
    else
        % Caso Giro 3D (phidot > 0): r es el radio polar siempre positivo
        ode_fun = @(t, Y) [ ...
            Y(2); ...
            ( Y(1)*( (p_phi / (m * max(Y(1), 1e-5)^2))^2 ) - 4*a^2*Y(1)*(Y(2)^2) - 2*g*a*Y(1) ) / (1 + 4*a^2*Y(1)^2); ...
            p_phi / (m * max(Y(1), 1e-5)^2) ...
        ];
    end
    
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    [t_ode, Y_ode] = ode45(ode_fun, tspan, [r; vr; phi], options);
    
    % --- Interpolación a tiempo uniforme (60 FPS) ---
    fps = 60; 
    num_frames = t_final * fps;
    t = linspace(0, t_final, num_frames);
    
    r_sol   = interp1(t_ode, Y_ode(:, 1), t);
    vr_sol  = interp1(t_ode, Y_ode(:, 2), t);
    phi_sol = interp1(t_ode, Y_ode(:, 3), t);
    
    % --- Cálculo Aceleración Radial y Normal ---
    if phidot == 0
        d2r_sol = (- 4*a^2*r_sol.*(vr_sol.^2) - 2*g*a*r_sol) ./ (1 + 4*a^2*(r_sol.^2));
    else
        d2r_sol = ( r_sol .* ((p_phi ./ (m * max(r_sol, 1e-5).^2)).^2) - 4*a^2*r_sol.*(vr_sol.^2) - 2*g*a*r_sol ) ./ (1 + 4*a^2*(r_sol.^2));
    end
    
    z_sol = a * (r_sol.^2);
    d2z_sol = 2 * a * (vr_sol.^2 + r_sol .* d2r_sol);
    
    % Magnitud analítica de la Reacción Normal N
    N_mag = abs(m * (d2z_sol + g) .* sqrt(1 + 4*a^2*(r_sol.^2)));
    % Coordenadas Cartesianas de la Masa
    x_sol = r_sol .* cos(phi_sol);
    y_sol = r_sol .* sin(phi_sol);
    
    % --- Configuración de la Figura ---
    fig = figure('Color', 'w', 'Name', 'Simulación Raíl Parabólico', 'Position', [100, 100, 900, 700]);
    ax = axes('Parent', fig);
    hold(ax, 'on'); grid(ax, 'off'); box(ax, 'off');
    view(ax, [-37.5, 20]);
    axis equal;
    xlim([-3.5, 3.5]); ylim([-3.5, 3.5]); zlim([-0.5, 4.8]);
    title('Atracción de Feria: La Langranbólica', 'Interpreter', 'latex', 'FontSize', 14);
    rotate3d on;
    
    % --- Ejes Fijos ---
    max_axis = 3.2;
    quiver3(0, 0, 0, max_axis, 0, 0, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.12);
    text(max_axis + 0.2, 0, 0, '$x$', 'Interpreter', 'latex', 'FontSize', 14);
    quiver3(0, 0, 0, 0, max_axis, 0, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.12);
    text(0, max_axis + 0.2, 0, '$y$', 'Interpreter', 'latex', 'FontSize', 14);
    quiver3(0, 0, 0, 0, 0, 4.2, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.12);
    text(0, 0, 4.5, '$z$', 'Interpreter', 'latex', 'FontSize', 14);
    text(-0.2, -0.2, -0.2, '$O$', 'Interpreter', 'latex', 'FontSize', 14, 'FontAngle', 'italic');
    
    % Raíl, Trayectoria y Masa
    r_rail_grid = linspace(-2.5, 2.5, 200);
    z_rail_grid = a * (r_rail_grid.^2);
    
    h_rail = plot3(r_rail_grid * cos(phi_sol(1)), r_rail_grid * sin(phi_sol(1)), z_rail_grid, 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5);
    h_trail = plot3(x_sol(1), y_sol(1), z_sol(1), 'r:', 'LineWidth', 1);
    h_particle = plot3(x_sol(1), y_sol(1), z_sol(1), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'k');
    text_m = text(x_sol(1) + 0.15, y_sol(1) + 0.15, z_sol(1), '$m$', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Vector Peso (mg)
    h_mg = quiver3(x_sol(1), y_sol(1), z_sol(1), 0, 0, -0.8, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
    text_mg = text(x_sol(1) + 0.1, y_sol(1), z_sol(1) - 0.8, '$m g$', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Vector Normal (N) - Verde (SOLO SI NO ROTA)
    if phidot == 0
        % Calculamos el vector inicial exacto para que quiver3 no se bugee
        norm_factor0 = sqrt(1 + 4*a^2*(r_sol(1)^2));
        nx0 = -2 * a * x_sol(1) / norm_factor0;
        ny0 = -2 * a * y_sol(1) / norm_factor0;
        nz0 = 1 / norm_factor0;
        scale0 = 0.0003 * N_mag(1) + 0.4;
        
        % Forzamos explícitamente 'AutoScale', 'off'
        h_N = quiver3(x_sol(1), y_sol(1), z_sol(1), nx0*scale0, ny0*scale0, nz0*scale0, ...
            'Color', 'g', 'LineWidth', 2.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
        text_N = text(x_sol(1) + 0.2, y_sol(1) + 0.2, z_sol(1) + 1.6, '$\vec{N}$', ...
            'Interpreter', 'latex', 'FontSize', 13, 'Color', [0 0.6 0]);
    end
    
    % --- Bucle de Animación ---
    dt_frame = 1 / fps;
    tic;
    for k = 1:length(t)
        phik = phi_sol(k);
        xk = x_sol(k);
        yk = y_sol(k);
        zk = z_sol(k);
        rk = r_sol(k);
        
        % 1. Actualizar rotación del raíl
        set(h_rail, 'XData', r_rail_grid * cos(phik), 'YData', r_rail_grid * sin(phik), 'ZData', z_rail_grid);
        
        % 2. Actualizar masa y estela
        set(h_trail, 'XData', x_sol(1:k), 'YData', y_sol(1:k), 'ZData', z_sol(1:k));
        set(h_particle, 'XData', xk, 'YData', yk, 'ZData', zk);
        set(text_m, 'Position', [xk + 0.15, yk + 0.15, zk]);
        
        % 3. Actualizar Peso mg
        set(h_mg, 'XData', xk, 'YData', yk, 'ZData', zk);
        set(text_mg, 'Position', [xk + 0.1, yk, zk - 0.8]);
        
        % 4. Actualizar Vector Normal N (SOLO SI NO ROTA)
        if phidot == 0
            norm_factor = sqrt(1 + 4*a^2*(rk^2));
            if abs(rk) < 1e-4
                nx = 0; ny = 0; nz = 1;
            else
                nx = -2 * a * xk / norm_factor;
                ny = -2 * a * yk / norm_factor;
                nz = 1 / norm_factor;
            end
            
            scale_visual = 0.0003 * N_mag(k) + 0.4; 
            Nx_v = nx * scale_visual;
            Ny_v = ny * scale_visual;
            Nz_v = nz * scale_visual;
            
            set(h_N, 'XData', xk, 'YData', yk, 'ZData', zk, ...
                     'UData', Nx_v, 'VData', Ny_v, 'WData', Nz_v);
            set(text_N, 'Position', [xk + Nx_v + 0.1, yk + Ny_v + 0.1, zk + Nz_v + 0.1]);
        end
        
        drawnow;
        tiempo_transcurrido = toc;
        tiempo_esperado = k * dt_frame;
        if tiempo_esperado > tiempo_transcurrido
            pause(tiempo_esperado - tiempo_transcurrido);
        end
    end
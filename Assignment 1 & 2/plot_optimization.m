function [figure_results] = plot_optimization(x, u, J, time, export_video)

export_video = false;

figure_results = figure('Name', 'optimization');
tiles = tiledlayout(1, 2, 'TileSpacing', 'tight');
% title(tiles, 'Space Shuttle Simulation');
% subtitle(tiles, title_string);

%% Tiles setup

for time_idx = 1:size(x, 2)

    nexttile(tiles, 1);
    s = stackedplot(time, [x{time_idx} u{time_idx}]);
    
    title('States and control optimization')
    xlabel("Time [s]")

    s.GridVisible = 1;
    s.FontSize = 14;
    s.DisplayLabels = ["x_1[m]", "v[m/s]", "m[kg]", "alpha[°]", "T[N]"];

    nexttile(tiles, 2);
    plot(0:length(J{time_idx})-1, J{time_idx})
    grid on
    
    title('Cost function J(x, u)')
    xlabel('Iteration #')

    drawnow
end



%% Animation export
if export_video
    video_object = VideoWriter(['submitted\' title_string '.mp4'], 'MPEG-4');
    video_object.Quality = 100;
    video_object.FrameRate = 2;
    
    open(video_object);
    writeVideo(video_object, plot_figures);
    close(video_object);
end

end


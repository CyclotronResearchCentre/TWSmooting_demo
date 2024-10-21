function aj_smoothing_plot(name, gsP_signal, twsP_signal, f_twsP_signal, tosP_signal, f_tosP_signal, ph_1D, dim)
% Function to visualize the results of different smoothing techniques (Gaussian, TWS, TSPOON) in 1D or 2D.
%
% INPUT:
% gsP_signal: Gaussian smoothed signal
% twsP_signal: Tissue-weighted smoothed signal (TWS)
% f_twsP_signal: Final tissue-weighted smoothed signal after applying mask
% tosP_signal: TSPOON smoothed signal per tissue
% f_tosP_signal: Final TSPOON smoothed signal after applying mask
% ph_1D: Original 1D signal
% dim: Data dimension (1D or 2D)
%
% OUTPUT:
% None (Plots are generated)

    switch dim
        case 1  % 1D Smoothing
            figure('Name', sprintf('1D results for file: %s', name));
            
            % Gaussian Smoothing
            subplot(3, 1, 1);
            plot(ph_1D, 'b-', 'DisplayName', 'Raw Data');
            hold on;
            plot(gsP_signal, 'r-', 'DisplayName', 'Gaussian Smoothing');
            %legend;
            xlabel('Position');
            ylabel('Signal Value');
            title('Gaussian Smoothing');

            % Tissue-Weighted Smoothing (TWS)
            subplot(3, 1, 2);
            plot(ph_1D, 'b-', 'DisplayName', 'Raw Data');
            hold on;
            plot(twsP_signal(1,:), 'r-', 'DisplayName', 'TWS (GM)');
            plot(twsP_signal(2,:), 'g-', 'DisplayName', 'TWS (WM)');
            plot(twsP_signal(3,:), 'y-', 'DisplayName', 'TWS (CSF)');
            plot(twsP_signal(4,:), 'Color', [1, 1, 0], 'DisplayName', 'TWS (Sculpt)');
            plot(f_twsP_signal, 'k-', 'DisplayName', 'TWS Final');
            %legend;
            xlabel('Position');
            ylabel('Signal Value');
            title('Tissue-Weighted Smoothing (TWS)');

            % TSPOON Smoothing
            subplot(3, 1, 3);
            plot(ph_1D, 'b-', 'DisplayName', 'Raw Data');
            hold on;
            plot(tosP_signal(1,:), 'r-', 'DisplayName', 'TSPOON (GM)');
            plot(tosP_signal(2,:), 'g-', 'DisplayName', 'TSPOON (WM)');
            plot(tosP_signal(3,:), 'm-', 'DisplayName', 'TSPOON (CSF)');
            plot(tosP_signal(4,:), 'Color', [1, 1, 0], 'DisplayName', 'TSPOON (Sculpt)');
            plot(f_tosP_signal, 'k-', 'DisplayName', 'TSPOON Final');
            %legend;
            xlabel('Position');
            ylabel('Signal Value');
            title('TSPOON Smoothing');

        case 2  % 2D Smoothing
            % Gaussian Smoothing
            figure('Name', sprintf('2D Gaussian Smoothing for file: %s', name));
            imagesc(gsP_signal);  % Display 2D data
            axis image;  % Maintain aspect ratio
            colorbar;  % Add color bar
            title('Gaussian Smoothing');
            xlabel('Position X');
            ylabel('Position Y');

            % TWS (Tissue-Weighted Smoothing) per tissue
            figure('Name', sprintf('2D Tissue-Weighted Smoothing for file: %s', name));
            
            % Display for GM tissue
            subplot(2, 2, 1);  % Create the first subplot in a 2x2 grid
            imagesc(squeeze(twsP_signal(1,:,:)));  % Display the smoothed signal for GM
            axis image;  % Ensure the aspect ratio is equal
            colorbar;  % Add a color bar
            title('Tissue-Weighted Smoothing (GM)');  % Title for GM tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Display for WM tissue
            subplot(2, 2, 2);  % Create the second subplot
            imagesc(squeeze(twsP_signal(2,:,:)));  % Display the smoothed signal for WM
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (WM)');  % Title for WM tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Display for CSF tissue
            subplot(2, 2, 3);  % Create the third subplot
            imagesc(squeeze(twsP_signal(3,:,:)));  % Display the smoothed signal for CSF
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (CSF)');  % Title for CSF tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Display for Sculpt tissue
            subplot(2, 2, 4);  % Create the fourth subplot
            imagesc(squeeze(twsP_signal(4,:,:)));  % Display the smoothed signal for Sculpt
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (Sculpt)');  % Title for Sculpt tissue
            xlabel('Position X');
            ylabel('Position Y');
            
            % Adjust the layout for better spacing
            sgtitle('Tissue-Weighted Smoothing');

            % TSPOON Smoothing
            figure('Name', sprintf('2D TSPOON Smoothing for file: %s', name));
            
            % Display for GM tissue
            subplot(2, 2, 1);  % Create the first subplot in a 2x2 grid
            imagesc(squeeze(tosP_signal(1, :, :)));  % Display the smoothed signal for GM
            axis image;  % Ensure the aspect ratio is equal
            colorbar;  % Add a color bar
            title('TSPOON Smoothing (GM)');  % Title for GM tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Display for WM tissue
            subplot(2, 2, 2);  % Create the second subplot
            imagesc(squeeze(tosP_signal(2, :, :)));  % Display the smoothed signal for WM
            axis image;
            colorbar;
            title('TSPOON Smoothing (WM)');  % Title for WM tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Display for CSF tissue
            subplot(2, 2, 3);  % Create the third subplot
            imagesc(squeeze(tosP_signal(3, :, :)));  % Display the smoothed signal for CSF
            axis image;
            colorbar;
            title('TSPOON Smoothing (CSF)');  % Title for CSF tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Display for Sculpt tissue
            subplot(2, 2, 4);  % Create the fourth subplot
            imagesc(squeeze(tosP_signal(4, :, :)));  % Display the smoothed signal for Sculpt
            axis image;
            colorbar;
            title('TSPOON Smoothing (Sculpt)');  % Title for Sculpt tissue
            xlabel('Position X');
            ylabel('Position Y');

            % Adjust the layout for better spacing
            sgtitle('Tissue-SPecific smOOthing compeNsated');
            
        case 3 % 3D Smoothing
            % Gaussian Smoothing
            figure('Name', sprintf('3D Gaussian Smoothing for file: %s', name));
            imagesc(gsP_signal(:,:,128/2));
            axis image;
            colorbar;
            title('Gaussian Smoothing');
            xlabel('Position X');
            ylabel('Position Y');
            
            % TWS (Tissue-Weighted Smoothing) per tissue
            figure('Name', sprintf('3D Tissue-Weighted Smoothing for file: %s', name));
            
            % Display for GM tissue
            subplot(2, 2, 1);
            imagesc(squeeze(twsP_signal(1,:,:,128/2)));
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (GM)');
            xlabel('Position X');
            ylabel('Position Y');

            % Display for WM tissue
            subplot(2, 2, 2);
            imagesc(squeeze(twsP_signal(2,:,:,128/2)));
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (WM)');
            xlabel('Position X');
            ylabel('Position Y');

            % Display for CSF tissue
            subplot(2, 2, 3);
            imagesc(squeeze(twsP_signal(3,:,:,128/2)));
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (CSF)');
            xlabel('Position X');
            ylabel('Position Y');

            % Display for Sculpt tissue
            subplot(2, 2, 4);
            imagesc(squeeze(twsP_signal(4,:,:,128/2)));
            axis image;
            colorbar;
            title('Tissue-Weighted Smoothing (Sculpt)');
            xlabel('Position X');
            ylabel('Position Y');
            
            % Adjust the layout for better spacing
            sgtitle('Tissue-Weighted Smoothing');

            % TSPOON Smoothing
            figure('Name', sprintf('3D TSPOON Smoothing for file: %s', name));
            
            % Display for GM tissue
            subplot(2, 2, 1);
            imagesc(squeeze(tosP_signal(1,:,:,128/2)));
            axis image;
            colorbar;
            title('TSPOON (GM)');
            xlabel('Position X');
            ylabel('Position Y');

            % Display for WM tissue
            subplot(2, 2, 2);
            imagesc(squeeze(tosP_signal(2,:,:,128/2)));
            axis image;
            colorbar;
            title('TSPOON (WM)');
            xlabel('Position X');
            ylabel('Position Y');

            % Display for CSF tissue
            subplot(2, 2, 3);
            imagesc(squeeze(tosP_signal(3,:,:,128/2)));
            axis image;
            colorbar;
            title('TSPOON (CSF)');
            xlabel('Position X');
            ylabel('Position Y');

            % Display for Sculpt tissue
            subplot(2, 2, 4);
            imagesc(squeeze(tosP_signal(4,:,:,128/2)));
            axis image;
            colorbar;
            title('TSPOON (Sculpt)');
            xlabel('Position X');
            ylabel('Position Y');

            % Adjust the layout for better spacing
            sgtitle('Tissue-SPecific smOOthing compeNsated');

        otherwise
            error('Data dimension is incorrect in aj_smoothing_plot');
    end
end

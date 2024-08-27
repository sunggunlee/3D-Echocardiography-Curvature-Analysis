clear;
% Read DICOM file
info = resampleDicom('p009npa.dcm');
data = info.data;
% Get dimensions
widthspan = info.widthspan;
width = info.width;
height = info.height;
depth = info.depth;
timevals = info.NumVolumes;
% Create main figure
fig = figure('Position', [100, 100, 1200, 800]); % Enlarge the figure
% Define central slices indices for each dimension
central_width = round(width / 2);
central_height = round(height / 2);
central_depth = round(depth / 2);
% Pre-allocate array to store LV volumes over time
lv_volumes = zeros(info.NumVolumes, 1);
% Loop through each time point
% Conversion factor from pixels to cm
pixel_to_cm = widthspan / width;
% Initialize arrays to store curvature and direction data
curvature_data = zeros(2, 8, 25); % Dimensions: curvatures per slice, slices, timepoints
error_data = strings(2, 8, 25); % Dimensions: curvature errors per slice, slices, timepoints
for t = 1:25
  
   % Display central slice
   subplot(3,2,1);
   imshow(squeeze(data(:, central_height, :, t))', []);
   title(['Center Height Slice: #' num2str(t)]);
  
   % Manually draw the long axis
   [x, z] = ginput(2); % Get user input   
   hold on;
   plot(x, z, 'bo', 'MarkerSize', 10);
   line(x, z, 'Color', 'b', 'LineWidth', 2);
   midpoint_x = ((x(1)+x(2))/2);
   midpoint_z = ((z(1)+z(2))/2);
   plot(midpoint_x, midpoint_z, 'go', 'MarkerSize', 20);
   hold off;
   % Translate to make midpoint coincide with origin (y,x,z orientation)
   data_time_t = squeeze(data(:, :, :, t));
   translated_data = imtranslate(data_time_t, [0, central_width - midpoint_x, central_depth - midpoint_z]);
   % Define angle of rotation
   theta = (tan((x(2) - x(1)) / (z(2) - z(1))))*180/pi; % Angle in degrees
   % Rotate around the y-axis (neg theta because CCW)
   rotated_data = imrotate3(translated_data, -theta, [1, 0, 0], 'linear', 'crop');
   % Translate back to original position
   translated_back_data = imtranslate(rotated_data, [0, -(central_width - midpoint_x), -(central_depth - midpoint_z)]);
   % Display rotated image
   subplot(3,2,2);
   final_heightslice_data = squeeze(translated_back_data(:, central_height, :))';
   imshow(final_heightslice_data, [0 255]);
   title('Rotated Image - Center Height');
   % Plot the new line and midpoint on the rotated image
   translated_x = x - midpoint_x;
   translated_z = z - midpoint_z;
   rotated_x = translated_x * cosd(theta) - translated_z * sind(theta);
   rotated_z = translated_x * sind(theta) + translated_z * cosd(theta);
   final_x = rotated_x + midpoint_x;
   final_z = rotated_z + midpoint_z;
   final_midpoint_x = ((final_x(1)+final_x(2))/2);
   final_midpoint_z = ((final_z(1)+final_z(2))/2);
   hold on;
   plot(final_x, final_z, 'ro', 'MarkerSize', 5);
   line(final_x, final_z, 'Color', 'r', 'LineWidth', 1);
   plot(final_midpoint_x, final_midpoint_z, 'ro', 'MarkerSize', 5);
   hold off;
   % Show top short-axis slices
   subplot(3,2,3);
   topsliceZ = round(final_z(1)+20);
   imshow(squeeze(translated_back_data(:, :, topsliceZ, 1))', [0,255]);
   title('Top Short Axis: #1');
   [x2, y2] = ginput(1); % Get user input   
   hold on;
   plot(x2, y2, 'gx', 'MarkerSize', 5);
   hold off;
   % Show bottom short-axis slices
   subplot(3,2,4);
   bottomsliceZ = round(final_z(2)-20);
   imshow(squeeze(translated_back_data(:, :, bottomsliceZ, 1))', [0,255]);
   title('Bottom Short Axis: #1');
   [x3, y3] = ginput(1); % Get user input   
   hold on;
   plot(x3, y3, 'gx', 'MarkerSize', 5);
   hold off;
   % Get top and bottom coordinates
   topcoord = [x2, y2, topsliceZ];
   botcoord = [x3, y3, bottomsliceZ];
   % Display width sliced image
   subplot(3,2,5);
   widthslice_data = squeeze(translated_back_data(central_width, :, :, 1))';
   imshow(widthslice_data, [0,255]);
   title('Center Width Slice: #1');
   hold on;
   plot(y2, topsliceZ, 'wx--', 'MarkerSize', 5)
   plot(y3, bottomsliceZ, 'wx--', 'MarkerSize', 5)
   y_new = [y2, y3];
   z_new = [topsliceZ, bottomsliceZ];
   line(y_new, z_new, 'Color', 'b', 'LineWidth', 1);
   new_mid_y = ((y_new(1)+y_new(2))/2);
   new_mid_z = ((z_new(1)+z_new(2))/2);
   plot(new_mid_y, new_mid_z, 'wo', 'MarkerSize', 5);
   hold off;
   % Translate to make midpoint coincide with origin (y,x,z orientation)
   data2_time_t = squeeze(translated_back_data(:, :, :, 1));
   translated_data2 = imtranslate(data2_time_t, [central_height - new_mid_y, 0, central_depth - new_mid_z]);
   % Define angle of rotation
   theta2 = -(atan((y_new(2) - y_new(1)) / (z_new(2) - z_new(1))))*180/pi; % Angle in degrees
   % Rotate around the y-axis (neg theta because CCW)
   rotated_data2 = imrotate3(translated_data2, -theta2, [0, 1, 0], 'linear', 'crop');
   % Translate back to original position
   translated_back_data2 = imtranslate(rotated_data2, [-(central_height - new_mid_y), 0, -(central_depth - new_mid_z)]);
   % Display rotated image
   subplot(3,2,6);
   final_widthslice_data = squeeze(translated_back_data2(central_width, :, :))';
   final_sideslice_data = squeeze(translated_back_data2(:, central_height, :))';
   imshow(final_widthslice_data, [0 255]);
   title('Rotated Image - Center Width');
   % Loop through each of the 8 slices
   pause(2); % Pause for 4 seconds
   clf; % Clear the current figure
   close(gcf);
  
   figure;
   imshow(final_sideslice_data, [0 255], 'InitialMagnification', 'fit');
   title('Final Reoriented - Center Width');
   [vent_w, vent_h] = ginput(2);
   hold on;
   plot(vent_w, vent_h, 'bo', 'MarkerSize', 10);
   line(vent_w, vent_h, 'Color', 'b', 'LineWidth', 2);
   hold off;
   clf; % Clear the current figure
   close(gcf);
   % Initialize variables to store radii and volume
   allsliceZ = zeros(1,8);
  
   for slice = 1:8
       ventriclepad = 5;
       slice_height = ((vent_h(2)-ventriclepad) - (vent_h(1)+ventriclepad)) / 8;
       % Determine slice location
       sliceZ = round((vent_h(1)+ventriclepad*2) + (slice_height * (slice-1)));
       if slice == 1 % Use '==' for comparison
           for val = 1:8
               slice_h = ((vent_h(2)-ventriclepad) - (vent_h(1)+ventriclepad/2)) / 8;
               sliceZ_heights = round((vent_h(1)+ventriclepad*2) + (slice_height * (val-1)));
               allsliceZ(val) = sliceZ_heights;
           end
           figure;
           imshow(final_sideslice_data, [0 255], 'InitialMagnification', 'fit');
           title('Final Reoriented - Center Width');
          
           for sliceval = 1:8
               slice_yval = allsliceZ(sliceval);
               hold on;
               plot([vent_w(2)-30, vent_w(2)+30], [slice_yval, slice_yval], 'r-');
               hold off;
               hold on;
               text(vent_w(2)+50, slice_yval, ['slice #', num2str(sliceval)], 'FontSize', 8, 'Color', 'red');
               hold off;
           end
           pause(3); % Pause for 4 seconds
           clf; % Clear the current figure
           close(gcf);
       end
       % Display short-axis slice
       short_axis_slice = squeeze(translated_back_data2(:, :, sliceZ, 1))';
       figure;
       imshow(short_axis_slice, [0, 255], 'InitialMagnification', 'fit');
       title(['Short Axis Slice: #' num2str(slice) '/8']);
     
       hold on;
       [curvature1, direction1, error1] = drawfreehand_curvature_direction1();
       hold off;
       hold on;
       [curvature2, direction2, error2] = drawfreehand_curvature_direction2(); % Get second curvature and direction
       hold off;
       % Store curvature and direction data for both readings
       curvature_data(1, slice, t) = curvature1;
       curvature_data(2, slice, t) = curvature2;
       error_data(1, slice, t) = error1;
       error_data(2, slice, t) = error2;
      
       pause(2); % Pause for 2 seconds
       clf; % Clear the current figure
       close(gcf);
      
   end
end
function [k_cm, direction, err] = drawfreehand_curvature_direction1()
   % Draw curve and get x,y coordinates
   hFH = drawfreehand();
   xy = hFH.Position;
   delete(hFH);
   hold on;
   xCoordinates = xy(:, 1);
   yCoordinates = xy(:, 2);
   plot(xCoordinates, yCoordinates, 'r.', 'LineWidth', 1, 'MarkerSize', 5);
   hold off;
   % Check inputs
   x = xy(:, 1);
   y = xy(:, 2);
   lx = length(x);
   if ~isvector(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x)) || ...
      ~isvector(y) || ~isfloat(y) || ~isreal(y) || ~all(isfinite(y)) || ...
       lx ~= length(y) || lx < 3
       error('Invalid input data.');
   end
   % Check collinearity
   if rank(diff([x([1:min(50,lx) 1]) y([1:min(50,lx) 1])])) == 1
       if lx <= 50 || rank(diff([x y;x(1) y(1)])) == 1
           k_cm = 0;
           err = NaN;
           direction = 'undefined';
           fprintf('Radius of curvature: %f\n', k_cm);
           fprintf('Root mean squared error: %f\n', err);
           fprintf('Direction of curvature: %s\n', direction);
           return;
       end
   end
   % Calculate curvature
   x = x(:);
   y = y(:);
   xx = x .* x;
   yy = y .* y;
   xy = x .* y;
   xxyy = xx + yy;
   sx = sum(x);
   sy = sum(y);
   sxx = sum(xx);
   syy = sum(yy);
   sxy = sum(xy);
   [L,U]=lu([sx sy lx;sxy syy sy;sxx sxy sx]);
   a=U\(L\[sxx+syy;sum(xxyy.*y);sum(xxyy.*x)]);
   xc=0.5*a(1);                % X-position of center of fitted circle
   yc=0.5*a(2);                % Y-position of center of fitted circle
   k=1/sqrt(xc^2+yc^2+a(3));	% Curvature of fitted circle
   k_cm = k * (299/18.8907);   % Convert curvature to centimeters
   % RMSE
   err = sqrt(mean((1./sqrt((x-xc).^2+(y-yc).^2)-k).^2));
   % Circle radius
   circle_radius = 1 / k;
   % Determine direction of curvature
   if x(4) < x(1)
       direction = -1; % Second point is to the left of the first point
   else
       direction = 1; % Second point is to the right of the first point
   end
   % Overlay circle with calculated radius of curvature
   hold on;
   theta = linspace(0, 2*pi, 100);
   x_circle = xc + circle_radius * cos(theta);
   y_circle = yc + circle_radius * sin(theta);
   plot(x_circle, y_circle, 'LineWidth', 0.5, LineStyle='--', Color=[1,0,0,0.25]);
   % Calculate midpoint of the line connecting first and last points
   midpoint_x = (x(1) + x(end)) / 2;
   midpoint_y = (y(1) + y(end)) / 2;
   % Calculate perpendicular line
   delta_x = x(end) - x(1);
   delta_y = y(end) - y(1);
   perp_x = midpoint_x + delta_y;
   perp_y = midpoint_y - delta_x;
   plot([midpoint_x, perp_x], [midpoint_y, perp_y], 'LineWidth', 0.25, Color=[0,0,1,0.25], LineStyle='--');
   % Plot lines with the same slope as the perpendicular line
   slope = -delta_x / delta_y;
   line_length = max(abs(x(end) - x(1)), abs(y(end) - y(1)));
   plot([x(1), x(1) + line_length], [y(1), y(1) + line_length * slope], 'LineWidth', 0.5, Color=[1,0,0,0.5], LineStyle='--');
   plot([x(end), x(end) + line_length], [y(end), y(end) + line_length * slope], 'LineWidth', 0.5, Color=[1,0,0,0.5], LineStyle='--');
   hold off;
   % Display results
   if direction == -1
       fprintf('Direction of curvature: Left\n');
   elseif direction == 1
       k_cm = -k_cm; % Convert curvature to negative for left curvature
       fprintf('Direction of curvature: Right\n');
   end
   fprintf('Radius of curvature: %f\n', k_cm);
   fprintf('Root mean squared error: %f\n', err);
end
function [k_cm, direction, err] = drawfreehand_curvature_direction2()
   % Draw curve and get x,y coordinates
   hFH = drawfreehand();
   xy = hFH.Position;
   delete(hFH);
   hold on;
   xCoordinates = xy(:, 1);
   yCoordinates = xy(:, 2);
   plot(xCoordinates, yCoordinates, 'b.', 'LineWidth', 1, 'MarkerSize', 5);
   hold off;
   % Check inputs
   x = xy(:, 1);
   y = xy(:, 2);
   lx = length(x);
   if ~isvector(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x)) || ...
      ~isvector(y) || ~isfloat(y) || ~isreal(y) || ~all(isfinite(y)) || ...
       lx ~= length(y) || lx < 3
       error('Invalid input data.');
   end
   % Check collinearity
   if rank(diff([x([1:min(50,lx) 1]) y([1:min(50,lx) 1])])) == 1
       if lx <= 50 || rank(diff([x y;x(1) y(1)])) == 1
           k_cm = 0;
           err = NaN;
           direction = 'undefined';
           fprintf('Radius of curvature: %f\n', k_cm);
           fprintf('Root mean squared error: %f\n', err);
           fprintf('Direction of curvature: %s\n', direction);
           return;
       end
   end
   % Calculate curvature
   x = x(:);
   y = y(:);
   xx = x .* x;
   yy = y .* y;
   xy = x .* y;
   xxyy = xx + yy;
   sx = sum(x);
   sy = sum(y);
   sxx = sum(xx);
   syy = sum(yy);
   sxy = sum(xy);
   [L,U]=lu([sx sy lx;sxy syy sy;sxx sxy sx]);
   a=U\(L\[sxx+syy;sum(xxyy.*y);sum(xxyy.*x)]);
   xc=0.5*a(1);                % X-position of center of fitted circle
   yc=0.5*a(2);                % Y-position of center of fitted circle
   k=1/sqrt(xc^2+yc^2+a(3));	% Curvature of fitted circle
   k_cm = k * (299/18.8907);   % Convert curvature to centimeters
   % RMSE
   err = sqrt(mean((1./sqrt((x-xc).^2+(y-yc).^2)-k).^2));
   % Circle radius
   circle_radius = 1 / k;
   % Determine direction of curvature
   if x(4) < x(1)
       direction = -1; % Second point is to the left of the first point
   else
       direction = 1; % Second point is to the right of the first point
   end
   % Overlay circle with calculated radius of curvature
   hold on;
   theta = linspace(0, 2*pi, 100);
   x_circle = xc + circle_radius * cos(theta);
   y_circle = yc + circle_radius * sin(theta);
   plot(x_circle, y_circle, 'LineWidth', 0.5, LineStyle='--', Color=[0,0,1,0.25]);
   % Calculate midpoint of the line connecting first and last points
   midpoint_x = (x(1) + x(end)) / 2;
   midpoint_y = (y(1) + y(end)) / 2;
   % Calculate perpendicular line
   delta_x = x(end) - x(1);
   delta_y = y(end) - y(1);
   perp_x = midpoint_x + delta_y;
   perp_y = midpoint_y - delta_x;
   %plot([midpoint_x, perp_x], [midpoint_y, perp_y], 'LineWidth', 0.25, Color=[0,0,1,0.25], LineStyle='--');
   % Plot lines with the same slope as the perpendicular line
   slope = -delta_x / delta_y;
   line_length = max(abs(x(end) - x(1)), abs(y(end) - y(1)));
   %plot([x(1), x(1) + line_length], [y(1), y(1) + line_length * slope], 'LineWidth', 0.25, Color=[1,0,0,0.25], LineStyle='--');
   %plot([x(end), x(end) + line_length], [y(end), y(end) + line_length * slope], 'LineWidth', 0.5, Color=[1,0,0,0.5], LineStyle='--');
   hold off;
   % Display results
   if direction == -1
       fprintf('Direction of curvature: Left\n');
       k_cm = -k_cm; % Convert curvature to negative for left curvature
   elseif direction == 1
       fprintf('Direction of curvature: Right\n');
   end
   fprintf('Radius of curvature: %f\n', k_cm);
   % fprintf('Root mean squared error: %f\n', err);
end

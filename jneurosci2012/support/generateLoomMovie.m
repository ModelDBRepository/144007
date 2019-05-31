function [M, t, theta, vel] = generateLoomMovie(t, loverv, degperpx, calibration_vect)
%function M = generateLoomMovie(t, loverv, degperpx)
% Generate a looming stimulus, the actual 2D stimulus movie
% t is the time vector, given in ms (time before collision is negative, approaching 0, collision-time)
% loverv - looming parameter, also in ms (should be positive)
% degperpx - The total size of the image is 180 deg in each direction, and this
%           parameter determines the size of each frame of the movie.
% The looming stimulus returned is black on white, and contains no
% anti-aliasing.
% The matrix returned is 4D, X x Y x 1 x T

%loverv = -40; %ms
%t = -1000:5:-1;
max_theta = 75;
theta = atan(loverv./t)*180/pi; %this is the half angle.
vel = -loverv./(t.^2+loverv^2)*180/pi;
max_angle_i = find(theta >= max_theta,1,'first');
if(~isempty(max_angle_i))
    theta(max_angle_i:end) = max_theta;
    vel(max_angle_i:end) = 0;
end

% slightly more sophisticated way to do the looming stimulus 
dt = mean(diff(t))/1000;
half_size = 400; %mm, simulated object physical half size
stim = loom2(dt, half_size/abs(loverv)); %generates the looming size at each frame.
t = -stim.t(end:-1:1) * 1000; %put the time into ms, reverse the time vector because the returned vector is backwards
theta = stim.loom_deg(end:-1:1, 2); %same here, reverse vector
centerx = 90/degperpx;
centery = 90/degperpx;
max_x = ceil(180/degperpx);
max_y = ceil(180/degperpx);
min_x = 0;
min_y = 0;
x = 0:max_x;
y = 0:max_y;

cmap = colormap(gray);
M = zeros(length(x), length(y), 1, length(t), 'single');
%M = single(M);
for i = 1:length(t)
    if (t(i) < 0)
        left = ceil(centerx - theta(i)/degperpx);
        lrem = (centerx-theta(i)/degperpx) - left; %use remainders for anti-aliasing
        left = max(2, left);

        right = floor(centerx + theta(i)/degperpx);
        rrem = (centerx + theta(i)/degperpx) - right;
        right = min(right, max_x-1);

        top = ceil(centery - theta(i)/degperpx);
        trem = (centery - theta(i)/degperpx) - top;
        top = max(2, top);

        bottom = floor(centery + theta(i)/degperpx);
        brem = (centery + theta(i)/degperpx) - bottom;
        bottom = min(bottom, max_y-1);  

        poly_x = [left-.1 left-.1 right+.1 right+.1]; %the fractions make sure that pixels on edge are included 
        poly_y = [bottom+.1 top-.1 top-.1 bottom+.1];
        loom_mask = poly2mask(poly_x, poly_y, length(x), length(y));
        neg_mask = 1-loom_mask;
        % need to put in some anti-aliasing, because the temporal resolution
        % at single pixels is really coarse without it. When indexing directly into the mask, 
        % vertical position is given by the first coordinate (since it is the row).
        neg_mask(top-1, left-1:right+1) = 1 - abs(trem);
        neg_mask(bottom+1,left-1:right+1) = 1- abs(brem);
        neg_mask(top-1:bottom+1, left-1) = 1- abs(lrem);
        neg_mask( top-1:bottom+1,right+1) = 1- abs(rrem);
        loom_im = neg_mask;
        % adjust the values to correspond to 1-256 values
        curr_frame = floor(63*loom_im + 1);
        % range checking
        bi = find(curr_frame < 1);
        if (~isempty(bi)) curr_frame(bi) = 1; end
        bi = find(curr_frame > 64);
        if (~isempty(bi)) curr_frame(bi) = 64; end
        calibframe = calibration_vect(curr_frame);
        M(:,:,1,i) = calibframe; %put in image sequence form  
    else
        M(:,:,1,i) = calibration_vect(1);
        
    end
   
end
M = squeeze(M(:,:,1,:));


% -----------------------------------------------------
function stim = loom2(Frame_Time,v, half_square_size)
% stim = loom2(Frame_Time,v)
% returns a stimulus structure with parameters of the stimulus
% and data to generate it
%
% halfsquaresize = size of the square from center to edge in mm, default is 400 mm
% v = velocity of approach in m/s
% halfmaxsize = maximal size from center to edge of the square in deg
% scr_an = distance screen animal (optional parameter) in mm
% output data: loom_deg and loom_pix
%
% Computes the size of a square recessing away from the animal
% at a constant velocity  on the screen from its
% maximal size down to a size of less than one deg. in visual angle.
% The size of a square is decreased from one frame to the next
% in such a way as to represent the projection on the screen of a square
% moving with constant velocity in the direction normal to the screen
% plane. The movement is recession only, the parameters are set
% in such a way that we can give real physical parameters. The corresponding
% approaching sequence can be obtained by performing a reverse sort (sort -n) to give an
% increasing sequence. The maximal size of the object on the screen (in
% degrees) is  given as a parameter and the recession proceeds until the square
% has an angular subtense of less than one degree.
% 
% Parameters: halfsquaresize v Frame_Time halfmaxdeg scr_an
%
% with halfsquaresize: The size of the square from its center to the edge in cm.
%                      This is effectively the half size of the square.
%
%      v:              Velocity of approach of the square in m/s.
%
%      Frame_Time:     Refresh time in ms
%            
%      halfmaxdeg:     Maximal size subtended by the object from the center to
%                      its edge on the screen in degree. This is effectively
%                      half the angular size of the object in degrees.
%
%      scr_an:         Distance between the eye and the animal in cm
%
%
% Remark: the last parameter scr_an is optional. If it is not given,
%         it is replaced by the standard value Screen_Animal defined below.
%         All the values are effectively computed for the half object,
%         since this turns out to be the easiest way to implement the
%         object motion on the screen.
%
% The program outputs two arrays:
%
% loom_deg = (t, y) coordinates of time (in sec) prior to collision
%            and half-angle of the object (angle center-edge)
%            in deg.
% loom_pix = (t, y) coordinates of time (in sec) prior to collision
%            and half angle of the object (angle center-edge)
%            in pixels
%
% l_v is the ratio of the square size divided by the speed which uniquely
% determines the sequence of images.
%
% These two files also contain as a comment the distance screen animal and
% the angular parameter determining the sequence of sizes.
%
% On the one hand, it is convenient to output the half-angle in degrees since
% this data can then be used for theoretical calculations and to determine
% the angular increase of an edge between two frames. On the other hand, the
% size in pixels is the data that is needed to program the visual stimulus.
% 

% 09/04/2008: Modified from lgeloom2.c in /bcm/gabbiani/lgmd 
% 09/09/2008: modified to test the accuracy of the results with 
%             files copied from saturn:lgmd/lge_stimulation/serie1_stims
% 06/11/2009: Takes frame rate as an argument

% Hardware values for the screen specific to the LGE system
% 
Horizontal_Pixel_Resolution = 640; %pixels on the screen
Horizontal_Size = 390; %corresponding length in mm
Vertical_Pixel_Resolution = 480;
Vertical_Size = 290;
x_as = 158; %distance screen animal in mm
%maximal size from center to edge of the square in deg, set to full screen
half_ang_max_h = floor(atand(0.5*Horizontal_Size/x_as));
half_ang_max_v = floor(atand(0.5*Vertical_Size/x_as));
half_ang_max = min(half_ang_max_h, half_ang_max_v);
% The use of vertical and not horizontal aspect ratio is arbitrary
px_per_mm = Vertical_Pixel_Resolution/Vertical_Size;
px_per_deg = .5 * Vertical_Pixel_Resolution/atand(0.5*Vertical_Size/x_as); %pixels/deg
mm_per_deg = Vertical_Size/atand(0.5*Vertical_Size/x_as);

if (~exist('half_square_size')) half_square_size = 400; end %size of the square from center to edge in mm, 400
v_mm_s = 1000*v; %in mm/s

%half_ang_max = floor(atan(min(Horizontal_Size, Vertical_Size)/2/x_as))

% SET KINEMATIC VALUES
dt = Frame_Time;  %frame time in seconds
d = half_square_size; %half square size in mm
l_v = d/v; %kinematic parameter in mm/(mm/ms) = ms 

% half the size of the object on the screen at the end of approach in mm
y_fin = x_as*tand(half_ang_max);

% COMPUTE IMAGE SEQUENCE
% 
% final distance of the object to the animal in mm
x_fin = d/tand(half_ang_max);
 
% final time remaining to collision
t_fin = x_fin/v_mm_s; %in sec
 
x = x_fin; %start at the final distance to the eye in mm
t = t_fin; %start at the final time to collision in sec
half_ang_size = half_ang_max; %final half angular size in deg
y = y_fin; %final half size of the object on the screen size in mm
px = round(y*px_per_mm); %half size of the object on the screen size in pixels

loom_deg(1,1:2) = [-t, half_ang_size];
loom_px(1,1:2) = [-t, px];
loom_px2(1, 1:2) = [-t, px];
loom_mm(1,1:2) = [-t, y];
loom_mm2(1,1:2) = [-t, y];
loom_x(1) = x;
t_vect(1) = t;
count = 1; 
 
while ( half_ang_size > 0.5 )
    count = count + 1;
    t = t + dt;
    x = x + dt*v_mm_s;
    half_ang_size = atand(d/x);
    y = (d/x)*x_as;
    px = round(y*px_per_mm);
    loom_deg(count,1:2) = [-t, half_ang_size];
    loom_px(count,1:2) = [-t, px];
    px2 = round(half_ang_size*px_per_deg);
    loom_px2(count,1:2) = [-t, px2];
    loom_mm(count, 1:2) = [-t, y];
    y2 = half_ang_size*mm_per_deg;
    loom_mm2(count,1:2) = [-t y2];
    loom_x(count) = x;
    t_vect(count) = t;
end;

n_frames = length(loom_px);
loom_sqr = zeros(n_frames,4);
loom_sqr(n_frames:-1:1,3) = 2*loom_px(:,2);
loom_sqr(n_frames:-1:1,4) = 2*loom_px(:,2);

stim.loom_deg = loom_deg;
stim.loom_px = loom_px;
stim.loom_px2 = loom_px2;
stim.loom_sqr = loom_sqr;
stim.t = t_vect; t - t_fin;
stim.l_v = l_v;
stim.loom_mm = loom_mm;
stim.loom_mm2 = loom_mm2;
stim.x = loom_x;

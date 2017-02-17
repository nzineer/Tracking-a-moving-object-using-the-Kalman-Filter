
           
%{

clear all;
close all;
clc;
%---------------------------------------------------Working on Centroid And
%Storing All of The Centroid values in a matrix 
cent=[];
vid=VideoReader('temp.mp4');
d=read(vid);
[r,c,t,n]=size(d);
for o = 1:n
    im=d(:,:,:,o);
    total_pixel=0;
    for x=1:10
        bw=zeros(r,c);
        for k=1:r
            for l=1:c   
                if(im(k,l,3)-im(k,l,1)>70 && im(k,l,3)-im(k,l,2)>70)
                    bw(k,l)=1;
                    total_pixel=total_pixel+1;
                end
            end
        end
    end
    imshow(bw);
    impixelinfo;
    s=strel('disk',3);
    bw=imerode(bw,s);
    bw = bwareaopen(bw, 200);
    centroid  = regionprops(bw,'Centroid');
    %v=size(centroid.Centroid);
    if total_pixel<=5000
        disp('Not found');
        cent=[cent;-1 -1];
    end
    cent = [cent; centroid.Centroid(1) centroid.Centroid(2)]; 
end
%----------------------------------Centroids have been
%Extracted+Stored+Printed
%}
%% define main variables
data=[];
dt = 1;  %our sampling rate
S_frame = 10; %starting frame
u = .005; % define acceleration magnitude
Q= [cent(1,1); cent(1,2); 0; 0]; %initized state--it has four components: [positionX; positionY; velocityX; velocityY] of the hexbug
Q_estimate = Q;  %estimate of initial location estimation of where the hexbug is (what we are updating)
Accel_noise_mag = 10; %process noise: the variability in how fast the Hexbug is speeding up (stdv of acceleration: meters/sec^2)
measurement_noise_x = 1;  %measurement noise in the horizontal direction (x axis).
measurement_noise_y = 1;  %measurement noise in the horizontal direction (y axis).
Ez = [measurement_noise_x 0; 0 measurement_noise_y];
Ex = [dt^4/4 0 dt^3/2 0; ...
    0 dt^4/4 0 dt^3/2; ...
    dt^3/2 0 dt^2 0; ...
    0 dt^3/2 0 dt^2].*Accel_noise_mag^2; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial Hexbug position variance (covariance matrix)

%% Define update equations in 2-D! (Coefficent matrices): A physics based model for where we expect the HEXBUG to be [state transition (state + velocity)] + [input control (acceleration)]
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrice
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0];  %this is our measurement function C, that we apply to the state estimate Q to get our expect next/new measurement


%% initize result variables
% Initialize for speed
Q_loc = []; % ACTUAL hexbug motion path
vel = []; % ACTUAL hexbug velocity
Q_loc_meas = []; % the hexbug path extracted by the tracking algo

%% initize estimation variables
Q_loc_estimate = []; %  position estimate
vel_estimate = []; % velocity estimate
P_estimate = P;
predic_state = [];
predic_var = [];
r = 5; % r is the radius of the plotting circle
j=0:.01:2*pi; %to make the plotting circle
vid=VideoReader('temp.mp4');
d=read(vid);

for t =1:n
    %% do the kalman filter  
    im=d(:,:,:,t);
    
    if(cent(t,1)~=-1 && cent(t,2)~=-1)
        Q_loc_meas(:,t) = [cent(t,1); cent(t,2)];
        
        Q_estimate(3)=-9;
        Q_estimate(4)=-9;
        u=30; 
    else
       % df=size(predic_state);
        Q_loc_meas(:,t)=[Q_estimate(1) Q_estimate(2)];
        Q_estimate(3)=-9;
        Q_estimate(4)=-9;
        u=30;
        
    end
    % Predict next state of the Hexbug with the last state and predicted motion.
    Q_estimate = A * Q_estimate + B * u;
    predic_state = [predic_state; Q_estimate(1)] ;
    %predict next covariance
    P = A * P * A' + Ex;
    predic_var = [predic_var; P] ;
    % predicted Ninja measurement covariance
    % Kalman Gain
    K = P*C'*inv(C*P*C'+Ez);
    % Update the state estimate.
    if ~isnan(Q_loc_meas(:,t))
        Q_estimate = Q_estimate + K * (Q_loc_meas(:,t) - C * Q_estimate);
    end
    % update covariance estimation.
    P =  (eye(4)-K*C)*P;
    
    %% Store data
    Q_loc_estimate = [Q_loc_estimate; Q_estimate(1:2)];
    vel_estimate = [vel_estimate; Q_estimate(3:4)];
   
    %% plot the images with the  tracking
    imshow(im);
    axis off
    colormap(gray);
    hold on;
    data=[data; Q_estimate(1) Q_estimate(2)];
    plot(r*sin(j)+Q_loc_meas(1,t),r*cos(j)+Q_loc_meas(2,t),'.g'); % the actual tracking
    plot(r*sin(j)+Q_estimate(1),r*cos(j)+Q_estimate(2),'.r'); % the kalman filtered tracking
    hold off
    pause(0.05)
end




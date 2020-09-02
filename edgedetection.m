clc;clear;close all;

%% Initialization
%Read image 
A=imread('test05.jpg'); %Read image
szA = size(A);

% convert 3D rgb images to grayscale
if ( length(szA)>2 )
    A = rgb2gray(A);
end
    
figure(1),subplot(1,2,1), 
imshow(A), title('Original Image'); 

%% Task 1: Gaussian Blur
gaussFilt = (1/159)*[2 4 5 4 2; 4 9 12 9 4; 5 12 15 12 5; 4 9 12 9 4; 2 4 5 4 2];

B = convfilt(A,gaussFilt);

%Image without Noise after Gaussian blur
B = uint8(B);
figure(1), subplot(1,2,2),
imshow(B), title('Blurred Image');

%% Task 2 
sobel_filter=[-1 0 1; -2 0 2; -1 0 1]; %For X-direction
%sobel_y=[-1 -2 -1; 0 0 0; 1 2 1]; %For y direction

Gx = convfilt(B,sobel_filter);       % x-dir
Gy = convfilt(B,sobel_filter');     % y-dir
    

figure(3), 
subplot(1,2,1), imshow(mat2gray(Gx)),title('Vertical Edges'); %Check Gradient in x
subplot(1,2,2), imshow(mat2gray(Gy)),title('Horizontal Edges'); %Check gradient in y
%% Task 3 : Gradient Magnitude
G = sqrt(Gx.^2 + Gy.^2);   % magnitude 

figure(5), imshow(mat2gray(G)), title('Edge Magnitude');
%% Task 4 & 5: Gradient Orientation
% Set
%equal compare for non-maxima supression
equal = 0;

theta = atand(Gy./Gx);         
theta = mod(theta+22.5,360);
theta = mod(floor(theta/45),8);   % convert to 8 cardinal dir
theta = mod(theta,4)+1;           % opposite directions are equal

thetacl = zeros([szA(1:2) 3]);    % init an rgb matrix for colour of the dir  
    
% Initialize matrix G for non-maximal suppression
G = round(G);  % round to nearest integer for comparison
Gnms = G;

% Make a padded G for use for non-maximal suppression
Gpad = padarray(G,[1,1],'symmetric');

% repeat over all pixels
for y=1:szA(1)
    for x=1:szA(2)

        direc = theta(y,x);  % direction

        if (~isnan(direc))

            % Prepares coloured illustration
            thetacl(y,x,:) = 255*de2bi(direc,3);    

            % Non-maximal suppression with thresholding
            xp = x+1;
            yp = y+1;

            switch (direc)
                case 1   % {0,180} degrees
                    cmp1 = Gpad(yp,xp-1);
                    cmp2 = Gpad(yp,xp+1);

                case 2   % {45,-135} degrees
                    cmp1 = Gpad(yp-1,xp-1);
                    cmp2 = Gpad(yp+1,xp+1);

                case 3   % {90,-90} degrees
                    cmp1 = Gpad(yp-1,xp);
                    cmp2 = Gpad(yp+1,xp);

                case 4   % {135,-45} degrees
                    cmp1 = Gpad(yp-1,xp+1);
                    cmp2 = Gpad(yp+1,xp-1);
            end



            comp = [Gpad(yp,xp), cmp1, cmp2];
            [~, maxIdx] = max(comp);

            if (equal)
                if (maxIdx ~= 1)
                    Gnms(y,x) = 0;
                end
            else
                if ~( (comp(1)>comp(2)) &&  (comp(1)>comp(3)) )
                    Gnms(y,x) = 0;
                end
            end
        end

    end
end

% Threshold and convert to binary
imTh = imbinarize(Gnms,multithresh(Gnms));
    
figure(6);
subplot(1,2,1), imshow(Gnms), title('Thinned edges');
subplot(1,2,2), imshow(imTh), title('Thinned edges with thresholding');

% Sets up the colourmap to show the gradient direction---------%
thval = linspace(0,315,8);  %
points = [cosd(thval); sind(thval)]';
p0 = [0 0];

figure(7), subplot(1,2,1);
for i=1:length(thval)
    p1 = points(i,:);
    dp = p1-p0;
    colour = de2bi(mod(i-1,4)+1,3);
    quiver(p0(1),p0(2),dp(1),dp(2),0, 'color', colour);
    hold on    
end
hold off
set(gca,'color','k');
set(gca, 'YDir','reverse')
title('Colourmap for edge orientation');

subplot(1,2,2),
imshow(thetacl), title('Gradient Orientation')


   
%% 1. Data Preparation

% Modes to include in summation:

W  = 10;        %Select The Womersley Cases To Analyse
Pa = 6000;          %Select Dimensionless Pressure Amplitude
n  = [1:1:20]; 
N  = length(n);    

%Time Discretization

dt  = 0.05; 
t   = [0:dt:15]; 
n_t = length(t); 

%Space Discretization

dy  = 0.05; 
y   = [-1:dy:1]; 
n_y = length(y);

%% 2. Inizialize Matrices :

Y   = zeros(n_y,N);      %Initialize Spatial Basis 
A_n = zeros(N,N);        %Initialize Amplitude Matrix 
T   = zeros(n_t ,N);     %Initialize Temporal Basis 
U_A = zeros(n_t,n_y);    %Initialize PDE solution

%% 1. Construct Spatial Basis

for i =1:1:length (n)
    N = 2 * n(i) - 1; % odd number in the series Y ( : , i ) = c o s ( N∗ p i ∗ y / 2 ) ;
end

% 2. Construct the Amplitudes

for i =1:1: length(n)
    N = 2 * n(i) - 1; % odd numbers in the series
    A_n(i ,i) = (16 * Pa) / (N * pi * sqrt((2 * W)^4 + N^4 * pi^4)); 
end

%% 3. Construct the Temporal Modes

for i = 1:1:length(n)
    N = 2 * n(i) - 1; % odd number in the series
    T(:,i) = (-1)^(n(i)) * cos( t - atan((4 * W^2)/(N^2 * pi^2))); 
end

%% 4. Assembly Solution

U_A = Y * A_n * T';

%% 5. Plot Result (ADDING BACK THE MEAN)

u_M = (1 - y.^2) * 0.5; % Compute the mean Flow u_M=repmat(u_Mb,length(t) ,1);%Repeat mean to obtain a matrix
u_T = U_A + u_M' ;       % Complete Solution

%% 3D PLOT ( With Function Print_3D_PLOT )

Fig1 = figure(1); 
surf(t,y,u_T)

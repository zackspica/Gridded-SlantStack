% SPICA DAS RECORDSECTION SCAN ALL EVENTS AND EXTRACT SW
% Spica et al., Subsurface Imaging with Ocean-Bottom Distributed Acoustic Sensing and Water Phases Reverberations
% Jorge Alberto Castillo Castellanos
% California Institute of Technology
% download readNPY.m @ https://github.com/kwikteam/npy-matlab   >> /Local_Tools/
% Disclamer : do the work, not optimized, take your time to run it. 


clear; close all; clc
addpath(genpath([pwd, '/Local_Tools/']))

%% Workpaths and input file
Workpath = 'Path/to/folder'; 
Datapath = [Workpath,'Data/'];

% Reference dispersion curve:
DispersionCurve = [Workpath, 'Files/Reference_DispersionCurve.txt'];

%% User input parameters
% Propagation direction:
Prop_Dir = 'Forward';                                                                                                                  % Forward or Backward

% Data grid parameters:
SR = 100;
Spacing = 5;
DIST_Vector = (2000 : Spacing : 10000)';
TIME_Vector = (0 : 1 / SR : 40)';

% Beamforming parameters:
Overlap = 90;                                                                                                                               % [percentage]
Vel_Per = 50;                                                                                                                                % [percentage]
VelocitySampling = 2;                                                                                                                 % [m/s]

%% Reading the reference dispersion curve
Curve = readtable(DispersionCurve, 'ReadVariableNames', 1, 'Delimiter', 'tab');
    
if ~strcmp(Curve.Properties.VariableNames{1}, 'Period')
    
    error('The reference dispersion curve must be in terms of period')
    
end

%% Defining the frequency limits of the band-pass filter cascade
Filter_Width = 0.07;                                                                                                                    % [Hz]
CenterFrequency = 1 ./ Curve.Period;                                                                                      % [Hz]
kfil = log10(CenterFrequency);                                                                                                 % [Hz]                          
fc1= 10.^(kfil - Filter_Width);                                                                                                    % [Hz]
fc2 = 10.^(kfil + Filter_Width);                                                                                                  % [Hz]
Filter_Corners = [fc1, fc2];                                                                                                         % [Hz] 

%% Listing all of the available directories
DIRList = dir(Datapath);

if isempty(DIRList)
    
    error('No .npz files found in the Datapath directory')
    
end

DIRList = extractfield(DIRList, 'name')';
DIRList = DIRList(~strncmpi('.', DIRList, 1));                                                                               %Remove every file starting with a "."

if isempty(DIRList)
    
    error('No files found in the Datapath directory')
    
end

%% Loading the data:
[XX, YY] = meshgrid(DIST_Vector, TIME_Vector);
Data_Matrix = nan(length(TIME_Vector), length(DIST_Vector), size(DIRList, 1));

for i = 1 : size(DIRList, 1)
    
    X = double(readNPY([Datapath, DIRList{i}, '/x.npy']));
    Y = double(readNPY([Datapath, DIRList{i}, '/y.npy']));
    Z = double(readNPY([Datapath, DIRList{i}, '/z.npy']));
    
    [X, Y] = meshgrid(X, Y);

    Data_Matrix(:,:,i) = interp2(X, Y, Z, XX, YY);
        
end

%% Looping through each event
[~, Max_Period_IDX] = max(Curve.Period);
DIST_WIN = Curve.Period(Max_Period_IDX) * Curve.Velocity(Max_Period_IDX);
TIME_WIN = 4 * DIST_WIN / Curve.Velocity(Max_Period_IDX);

% Overlapping time and distance vectors:
DIST_Overlap = DIST_WIN * ((100 - Overlap) / 100);
DIST_Windows = [(DIST_Vector(1) : DIST_Overlap : (DIST_Vector(end) - DIST_WIN))',...
                              ((DIST_Vector(1) + DIST_WIN) : DIST_Overlap : DIST_Vector(end))'];
DIST_Center = mean(DIST_Windows, 2);                    
TIME_Overlap = TIME_WIN * ((100 - Overlap) / 100);
TIME_Windows = [(TIME_Vector(1) : TIME_Overlap :...
                               (TIME_Vector(end) - TIME_WIN))',...
                              ((TIME_Vector(1) + TIME_WIN) :...
                              TIME_Overlap : TIME_Vector(end))'];
TIME_Center = mean(TIME_Windows, 2);

for i = 1 : size(DIRList, 1)
    
    Data = Data_Matrix(:,:,i);
    Data(isnan(Data)) = 0;
    
    % Initializing the image matrices:
    Prop_Velocity = nan(size(TIME_Center, 1), size(DIST_Center, 1), size(Filter_Corners, 1));
    Prop_Power = nan(size(TIME_Center, 1), size(DIST_Center, 1), size(Filter_Corners, 1));
    
    % Loopingth through the period ranges:
    for j = 1 : size(Filter_Corners, 1)
        
        % Reference period and velocity:
        Ref_Period = Curve.Period(j);
        Ref_Velocity = Curve.Velocity(j);

        % Velocity range:
        Min_Velocity = floor(Ref_Velocity - (Ref_Velocity * (Vel_Per / 100)));
        Max_Velocity = ceil(Ref_Velocity + (Ref_Velocity * (Vel_Per / 100)));
        Range_Velocity = (Min_Velocity : VelocitySampling : Max_Velocity)';
        
         % Filtering the data:
        Waveforms = zeros(size(Data));

        for k = 1 : size(Data, 2)

                Trace = Data(:,k);
                Trace_f = f_FiltSignal(Trace, 1 / SR, 8,...
                                                     Filter_Corners(j,1), Filter_Corners(j,2), 1, 'bandpass');       
                Waveforms(:,k) = Trace_f ./ max(abs(Trace_f));

        end
        
        for k = 1 : size(TIME_Center, 1)
        
            for kk = 1 : size(DIST_Center, 1)

                DIST_IDX = find(DIST_Vector >= DIST_Windows(kk,1) &...
                                            DIST_Vector <= DIST_Windows(kk,2));
                TIME_IDX = find(TIME_Vector >= TIME_Windows(k,1) &...
                                             TIME_Vector <= TIME_Windows(k,2));

                tPlotTimeSignals_Frame = TIME_Vector(TIME_IDX);
                Distance_Frame = DIST_Vector(DIST_IDX);
                Waveforms_Frame = Waveforms(TIME_IDX, :);
                Waveforms_Frame = Waveforms_Frame(:, DIST_IDX);
                
                if strcmp(Prop_Dir, 'Backward')
                    
                    Waveforms_Frame = fliplr(Waveforms_Frame);
                    
                end

                % Calculating the lag times:
                Initial_Time = tPlotTimeSignals_Frame(1);
                Initial_Offset = Distance_Frame(1);
                Distance_Frame = Distance_Frame - Initial_Offset;
                tau = zeros(size(Range_Velocity, 1), size(Distance_Frame, 1));

                for l = 1 : size(Distance_Frame, 1)

                    tau(:,l) = -Distance_Frame(l) ./ Range_Velocity;

                end

                tau_samples = round(tau ./ (1 / SR));

                % -------------------------------------------------------------------
                % Shift and Stack
                % -------------------------------------------------------------------
                Max_tau_Samples = max(abs(tau_samples(:)));
                Traces = zeros(size(Waveforms_Frame, 1) + (Max_tau_Samples * 2),...
                                           size(Waveforms_Frame, 2));
                Traces((Max_tau_Samples + 1) :...
                            (Max_tau_Samples + size(Waveforms_Frame, 1)),:) = Waveforms_Frame;

                 % Shifting the records:
                energy = zeros(size(tau_samples, 1), 1);

                for l = 1 : size(tau_samples, 1)

                    s = zeros(size(Traces));

                    for m = 1 : size(Traces, 2)

                        shift = tau_samples(l,m);
                        s(:,m) = circshift(Traces(:,m), shift);

                    end

                    energy(l) = sum(sum(s, 2).^2);

                end

                % Finding the best stacking velocity:
                [Max_energy, Max_IDX] = max(energy);
                Stacking_Velocity = Range_Velocity(Max_IDX);
                Predicted_Time = Distance_Frame ./ Stacking_Velocity;

                Prop_Velocity(k,kk,j) = Stacking_Velocity;
                Prop_Power(k,kk,j) = Max_energy;

            end
        
        end
        
        disp(['Frequency ', num2str(j), ' out of ', num2str(size(Filter_Corners, 1)),...
                 ' for event ', num2str(i), ' out of ', num2str(size(DIRList, 1)), ' processed'])
        
    end
    
    % Saving the output:
    Slant_Image.Periods = Curve.Period;
    Slant_Image.Velociy = Curve.Velocity;
    Slant_Image.Distance = DIST_Center;
    Slant_Image.Time = TIME_Center;
    Slant_Image.PowerImage = Prop_Power;
    Slant_Image.VelocityImage = Prop_Velocity;
    
    if strcmp(Prop_Dir, 'Forward')
        
        FileName = 'Slant_Image_Forward.mat';
        
    elseif strcmp(Prop_Dir, 'Backward')
        
        FileName = 'Slant_Image_Backward.mat';
        
    else
        
        error('Please select a forward or backwards propagation direction')
        
    end
    
    save([Datapath, DIRList{i}, '/', FileName], 'Slant_Image')
    
end
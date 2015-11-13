function [ phase_correction ] = probePhase(freqpts, use_case)

%experiment setup
distance=0.041; %Distance between the probe antennas for the probe response measurement
c=3*10^8; %Speed of light
measured=[];

switch use_case
    case 1
        %probe phase
        probe_Measurement_Path='D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Test Folder\Probe_Response_201_Freq.s2p';
        
        %Load the probe measurement in s2p format
        load_response=dlmread(probe_Measurement_Path,' ',9,0);
        
        load_response=load_response(:,[1 4 5]);
        probe_response=[load_response(:,1), load_response(:,2)+1.0j*load_response(:,3)];
        
        %choose desired frequency values
        [~, member]=ismember(freqpts,load_response(:,1));
        if length(load_response(:,1))<length(freqpts) || sum(ismember(member, 0))>0
            error('measurements do not cover frequency range of interest')
        end
        
        probe_response=probe_response(member,:);
        
        % Compensate for the distance between the probes
        measured=probe_response(:,2)./exp(-1.0j*distance*2*pi*probe_response(:,1)/c);
        
        % %compensate for connector used in calibrating cables before probe
        % %measurement. Assued to be cal kit through 85521A (115.881 ps delay)
        % delay_time=115.881e-12 %ps
        % measured=measured.*exp(-1.0j*2*pi*probe_response(:,1)*delay_time);
        
        % From two probes to one probe (Compensate for phase multiplicaction)
        phase_correction=abs(measured).*sqrt(exp(-1.0j*angle(measured)));
        
    case 2
        probe_Measurement_Path='C:\Users\lab\Documents\data\response\nsi\probephase.mat'; 
        load(probe_Measurement_Path)
        
        phase_correction=exp(1.0i*probe_phase_meas);
    case 3
        path='C:\Users\lab\Documents\data\response\nsi\nsi_wr42_18p0-26p0GHz_801pts.mat'; 
        probe_phase=load(path);
        
        [f_check, f_location]=ismember(freqpts,probe_phase.f);
        if ~all(f_check)
        error('Input frequencies do not match up with NSI data frequencies')
        end
        phase_correction=probe_phase.r(f_location);
        
end

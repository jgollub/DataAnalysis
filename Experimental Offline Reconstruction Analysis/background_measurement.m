


%%% initialize instraments
delete(instrfind) %delete any existing instrament objects 
vobj_switch = agilent_11713C_switchdriver_startVISAObject; %open switch communications
vobj_vna    = agilent_N5245A_NA_startVISAObject;           %open vna communications
calfile = 'SOLT_Cal_S31_S32_S21_p3Coupler_Off_6ft_13dBm_07082013';
[buffersize, F] = Initialize_N5245A_vers3(vobj_vna, 101, 17.5, 26.5,1E3,calfile); % setup VNA scan
pause(2) %wait for VNA 
num_switches=6; %six switches per bank
num_banks=2; %two banks (total of 12 switches)
for i_bank=1:num_banks
agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,repmat(i_bank,[1,num_switches]),[1:num_switches])
end

%% take bkg
% avg_factor = 50;
% %sweep through all switches
% background = zeros(length(F),num_switches,num_banks);
% for Switchbank=1:num_banks
%     switch Switchbank
%         case 1
%             sparam = 'MeasS31';
%         case 2
%             sparam = 'MeasS32';
%     end
%     
%     for sn1=1:num_switches
% 
%         agilent_11713C_switchdriver_openChannelNumbers(vobj_switch,Switchbank,sn1) %open switch
%         %measure s31 panelA to horn
%         for an=1:avg_factor
%             background(:,sn1,Switchbank) = background(:,sn1,Switchbank) + transpose(Read_N5245A(vobj_vna,buffersize,sparam));
%         end
%         background(:,sn1,Switchbank) = background(:,sn1,Switchbank)./avg_factor;
%         agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,Switchbank,sn1) %close switch
% 
%     end
% end

%% take empty chamber measurements
IFBW = [1E3,2E3,5E3, 10E3, 50E3, 100E3, 500E3, 1E6, 2E6, 3E6, 4E6, 5E6];
g_emptychamber = zeros(length(IFBW),length(F)*num_switches*num_banks);
instances = 50;

for nifbw=1:length(IFBW)
    ifbw = IFBW(nifbw)
    set_VNA_IFBandwidth(vobj_vna,ifbw)
    %sweep through all switches
    gempty = zeros(length(F),num_switches,num_banks,instances);
    
        %one full measurement vector
        for Switchbank=1:num_banks
        
            switch Switchbank
                case 1
                    sparam = 'MeasS31';
                case 2
                    sparam = 'MeasS32';
            end

            for sn1=1:num_switches
                agilent_11713C_switchdriver_openChannelNumbers(vobj_switch,Switchbank,sn1) %open switch
                for in=1:instances
                    %measure panel to horn
                    gempty(:,sn1,Switchbank,in) = transpose(Read_N5245A(vobj_vna,buffersize,sparam));
                end
                agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,Switchbank,sn1) %close switch
            end
        end
        
        for in=1:instances
            ginstance = gempty(:,:,:,in);
            g_emptychamber(nifbw,:,in) = transpose(ginstance(:)); 
        end
end


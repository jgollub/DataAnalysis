%{ 
                          TS April 2016
Comment block

Sample matlab program for the Keysight PNA vector network analyzers.
The sample program connects to a PNA Family network analyzer thru a VISA resource string. 

The program first clears the error queue and all status registers via the
"*CLS" command. The *IDN? identification query is then asserted and the resultant string is read.

The trigger is set to hold, after which a single trigger is asserted
with trigger completion hold-off via "INIT:CONT OFF" and "INIT:IMM;*OPC?"

Upon sensing of the trigger complete the binary bin block data transfer mode 
is set via "FORM:DATA REAL,64". For this 64 bit binary bin block transfer the byte order 
is swapped via the "FORM:BORD SWAP" command. 

Either the formatted data ("CALC:DATA? FDATA"), matching the display format, or alternatively 
the real and imaginary pairs of data("CALC:DATA? SDATA")representing the current trace data is
queried and returned (as a binary bin block, real 64 bit, data transfer). The data transfer 
method is then set back to ASCII via the "FORM:DATA ASCII".

As a wrap up the system error queue is checked at conclusion of the
application. If no errors were generated the response to the "SYST:ERR?"
query will still be "+0, "No Error"". 

%}



%Remove all interfaces to instrument
instrreset

% find all previously created objects
oldobjs = instrfind;

% If there are any existing objects
if (~isempty(oldobjs))
    % close the connection to the instrument
    fclose(oldobjs);
    % and free up the object resources
    delete(oldobjs);
end
 
% Remove the object list from the workspace.
clear oldobjs;

%Define PNA interface, this is the VISA resource string. Replace this VISA
%resource string with the your PNA VISA resource string as appropriate. 
pna = visa('agilent', 'TCPIP1::10.236.72.222::hpib7,16::INSTR');

%Buffer size must precede open command
set(pna,'InputBufferSize', 640000);
set(pna,'OutputBufferSize', 640000);

% open session to pna based on VISA resource string
fopen(pna);

% Clear the event status registers and all errors which may be in the PNA's queue.
fprintf(pna, '*CLS');

% Check to ensure the error queue is clear. Response is "+0, No Error"
fprintf(pna, 'SYST:ERR?'); 
errIdentifyStart = fscanf(pna, '%s');

%Query instrument identification string
fprintf(pna, '*IDN?'); 
idn = fscanf(pna, '%s');

%PNA timeout is set to 15 (seconds) to allow for longer sweep times. 
set(pna, 'Timeout', 15);

%Trigger mode is set to initiate continuous off
fprintf(pna, 'INIT:CONT OFF')

%Trigger a single sweep and wait for trigger completion via *OPC? (operation complete). 
fprintf(pna, 'INIT:IMM;*OPC?');
opComplete = fscanf(pna, '%s');

%Swap byte order on data query return.
fprintf(pna,'FORM:BORD SWAP');

%Set Trace Data read or return format as binary bin block real 64 bit values
fprintf(pna, 'FORM:DATA REAL,64');

%Select a trace to be read
fprintf(pna, 'CALC:PAR:SEL "CH1_S11_1"');

%To select the "FORMATTED DATA" which matches the display use the 
%'CALC:DATA? FDATA query. Alternatively, to select the underlying real and 
%imaginary pairs which the formatted data is based upon use 
%'CALC:DATA? SDATA query. Select or uncomment one of the following  
dataQueryType = 'CALC:DATA? FDATA' 
% dataQueryType = 'CALC:DATA? SDATA' 
fprintf(pna, dataQueryType);

%Read return data as binary bin block real 64 bit values. 
cData = binblockread(pna, 'float64')

% Binblock read has a 'hanging line feed that must be read and disposed
fscanf(pna, '%c')

%Read the stimulus values
fprintf(pna,'SENSE:X?');
frequency = binblockread(pna,'float64')
% Binblock read has a 'hanging line feed that must be read and disposed
fscanf(pna, '%c')

%Return data transfer format back to ASCII string format
fprintf(pna, 'FORM:DATA ASCII');

%As a last step requery the PNA error queue and ensure no errors have
%occured since intiation of program
fprintf(pna, 'SYST:ERR?'); 
errIdentifyStop = fscanf(pna, '%s');

%Close session connection
fclose(pna);
delete(pna)
clear pna

%Generate a plot of stimulus (frequency) - vs. - response (FDATA). 
if dataQueryType == 'CALC:DATA? FDATA'
        plot(frequency,cData) 
end 
%{
If the data is queuried as 'CALC:DATA? SDATA' then the returned array, 
as read in this application, is a 1xN array where N = 2 * Number of trace
points. Additional efforts, i.e. a routine, must be integrated to allow the 
resultant real and imaginary pairs to be parsed in to a 2 x N array. 
Then, the log magnitude is computed as 
    20*Log10*squareRoot(Real^2 + Imaginary^2). 
The X-Y, or Stimulus-Response can then be plotted viat the plot(x,y) 
invocation. 
%}
    

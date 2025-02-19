% set system limits
% Original gre.m modified by Didi
sys = mr.opts('MaxGrad',12,'GradUnit','mT/m','MaxSlew',100,'SlewUnit','mT/m/ms','rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6);

% Define FOV and resolution
fov = 256e-3;
sliceThickness = 5e-3;
Nx = 128;
Ny = Nx;

% Define sequence parameters
TE = 8e-3;
TR = 200e-3;
alpha=30;

% Create a new sequence object
seq=mr.Sequence(sys);

% Create slice selective alpha-pulse and corresponding gradients
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180, 'Duration', 4e-3,'SliceThickness', sliceThickness, 'apodization', 0.5,'timeBwProduct', 4, 'system' ,sys);
gz.amplitude = -gz.amplitude;
gzReph.amplitude = -gzReph.amplitude;

% Define other gradients and ADC events
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
% gx = mr.makeTrapezoid('x', 'FlatArea', -Nx*deltak, 'FlatTime', 6.4e-3,'system',sys);
gx = mr.makeTrapezoid('x', 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3,'system',sys);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', 2e-3,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% Calculate timing
delayTE = round((TE - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2 - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

spoilArea=4*gx.area(); % 4 "looks" good
% Add spoilers in read, refocus in phase and spoiler in slice
gxPost = mr.makeTrapezoid('x', 'Area', spoilArea, 'system', sys); % we pass 'system' here to calculate shortest time gradient
gyPost = mr.makeTrapezoid('y', 'Area', -max(phaseAreas(:)), 'Duration', 2e-3);
gzPost = mr.makeTrapezoid('z', 'Area', spoilArea, 'system', sys);

delayTR = delayTR - mr.calcDuration(gxPost, gyPost, gzPost);

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    seq.addBlock(rf, gz);
    gyPre = mr.makeTrapezoid('y', 'Area', phaseAreas(i), 'Duration', 2e-3);
    gyPre.amplitude = -gyPre.amplitude;
    gxPre.amplitude = -gxPre.amplitude;
    seq.addBlock(gxPre, gyPre, gzReph);
    seq.addBlock(mr.makeDelay(delayTE));
    gx.amplitude = -gx.amplitude;
    seq.addBlock(gx, adc);
    gyPost = mr.makeTrapezoid('y', 'Area', -gyPre.area, 'Duration', 2e-3);
    gyPost.amplitude = -gyPost.amplitude;
    % Add spoilers in read and slice and may be in phase
    gxPost.amplitude = -gxPost.amplitude;
    seq.addBlock(gxPost, gyPost, gzPost);
    seq.addBlock(mr.makeDelay(delayTR));
end

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
%%
% Trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold; plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.');


%%
% export definitions
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('TR', TR);
seq.setDefinition('TE', TE);
seq.setDefinition('BaseResolution',256);
seq.setDefinition('Name', 'DEMO_gre');

% seq.install('siemens');
% seq.write(['DEMO_gre_30_apod_p5.seq']);      % Write to pulseq file
seq.plot('timeRange', [0 2*TR]);

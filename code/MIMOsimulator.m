% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% MIMO Simulation Framework
% 1 User only (simplest case) 
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = MIMOsimulator(P)

    RX = P.RXperUser;
    TX = P.TXperUser;
    
    % Generate the spreading sequences
    HadamardMatrix = hadamard(P.HadLen);     
    
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence  = HadamardMatrix(:,randperm(P.HadLen,P.CDMAUsers)).';
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
               
    for user=1:P.CDMAUsers
        PNSequence(user,:) = step(LongCode).';
    end
    
    % Channel
    switch P.ChannelType
        case 'Multipath'
            NumberOfChipsRX = NumberOfChips+P.ChannelLength-1;
        otherwise
            NumberOfChipsRX = NumberOfChips;
    end

Results = zeros(1,length(P.SNRRange));

% Note that these are terminated, so they include the tail
trellis = poly2trellis(P.K, P.ConvSeq);
convEnc = comm.ConvolutionalEncoder(trellis, 'TerminationMethod', 'Terminated');
convDec = comm.ViterbiDecoder(trellis, 'TerminationMethod', 'Terminated', 'InputFormat','Hard');

for frame = 1:P.NumberOfFrames/P.CDMAUsers
    
    frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch P.Mode
        case 'HighRate'
            bits = randi([0 1],P.CDMAUsers,TX,P.NumberOfBits); % Random Data, High rate mode  
            % Add Frame Quality Indicator (bonus)
            bits_ind = cat(3,bits, randi([0 1],P.CDMAUsers,TX,P.Q_Ind));
        case 'HighDiversity'
            bits = randi([0 1],P.CDMAUsers,1,P.NumberOfBits); % Few data, redundancy
            % Duplicate for each antenna
            bits_dupl = ones(P.CDMAUsers,TX,1).*bits;
            % Add Frame Quality Indicator (bonus)
            bits_ind = cat(3,bits_dupl, randi([0 1],P.CDMAUsers,TX,P.Q_Ind));
        otherwise
            disp('Mode not supported');
    end

    % Convolutional encoding
    encoded_bits = zeros(P.CDMAUsers,TX,NbTXBits);
    for user=1:P.CDMAUsers
        for t=1:TX
            encoded_bits(user,t,:) = convEnc(squeeze(bits_ind(user,t,:))).';
        end
    end
    
    % Symbol repetition
    % Symbols shall not be repeated for a data rate of 9600bps, the one
    % we are simulating (7-6)
    % encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    
    % Interleaver
    if strcmp(P.Interleaving, 'On')
        for user=1:P.CDMAUsers
            encoded_bits(user,:,:)=matintrlv(squeeze(encoded_bits(user,:,:)).',P.IntrlvRows,NbTXBits/P.IntrlvRows).';
        end
    end
    
    % Pulse Shape (PNSequence)
    PN_symbols = zeros(P.CDMAUsers,TX, NbTXBits);
    for user = 1:P.CDMAUsers
        for t=1:TX
            PN_symbols(user,t,:) = xor(PNSequence(user,:), squeeze(encoded_bits(user,t,:)).'); % TODO add user loop here
        end
    end
    
    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1); 
    
    % Spreading with Hadamard
    symbol_spread = zeros(P.CDMAUsers,TX,P.HadLen, NbTXBits);
    for user=1:P.CDMAUsers
        for t = 1:TX
            symbol_spread(user,t,:,:) = HadamSequence(user,:) .* symbols(user,t,:);
        end
    end   
    waveform = reshape(symbol_spread,P.CDMAUsers,TX,NumberOfChips); 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = 1;
        case 'Multipath'
            himp = sqrt(1/2)* (randn(P.CDMAUsers,RX,TX,P.ChannelLength) +...
                            1i * randn(P.CDMAUsers,RX,TX,P.ChannelLength));
        otherwise
            disp('Channel not supported')
    end
    %Normalization with Euclidean norm along third axis
    himp = himp./vecnorm(himp, 2, 4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    
    switch P.ChannelType
        case {'ByPass', 'AWGN'}
            snoise = (randn(1,NumberOfChips) + ...
                      1i* randn(1,NumberOfChips) );
        case 'Multipath'
            snoise = (randn(P.CDMAUsers,RX,TX,NumberOfChipsRX) + ...
                      1i* randn(P.CDMAUsers,RX,TX,NumberOfChipsRX) );    
    end
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) * snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass' 
                y = waveform * himp;
            case 'AWGN'
                y = waveform * himp + noise;  
            case 'Multipath'
                y = zeros(P.CDMAUsers,RX,NumberOfChipsRX); %Normally add the users here!
                for userRX = 1:P.CDMAUsers
                    for r = 1:RX
                        for userTX = 1:P.CDMAUsers
                            for t = 1:TX
                                y(userRX,r,:) = squeeze(y(userRX,r,:)).' + conv(squeeze(waveform(userTX,t,:)),squeeze(himp(userTX,r,t,:))).' + squeeze(noise(userTX,r,t,:)).';
                            end
                        end
                    end
                end
            otherwise
                disp('Channel not supported')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Receiver
        switch P.ChannelType % Rake receiver makes no sense in this case
            case {'AWGN','ByPass'}
                % Despreading
                rxsymbols = conj(himp).*HadamSequence.'*reshape(y, P.HadLen, NbTXBits);
                desp_bits = rxsymbols < 0;
                
            case 'Multipath'
                % Despreading
                rxsymbols = zeros(P.CDMAUsers,RX,P.RakeFingers,NbTXBits); % TODO diversity or High rate possibility???
                himp_perm = zeros(P.CDMAUsers,RX,P.RakeFingers,TX);
                
                for user = 1:P.CDMAUsers
                    % Separation between antennas ? How to get himp ?
                    for r = 1:RX
                        % Order the best combination of fingers (all contribution
                        % are combiend and unseparable !
                        [~,ind] = maxk(sum(abs(himp(user,r,:,:)),3),P.RakeFingers,4);
                        ind = squeeze(ind);
                        % Permute corresponding elements of himp for later
                        % calculation
                        himp_perm(user,r,:,:) = permute(himp(user,r,:,ind), [1,2,4,3]);
                        
                        for finger = 1:P.RakeFingers
                            if strcmp(P.ChannelType,'Multipath')
                                % Despreading
                                rxsymbols(user,r,finger,:) = HadamSequence(user,:) * reshape(...
                                    y(user,r,ind(finger):ind(finger)+NumberOfChips-1),...
                                    P.HadLen, NbTXBits);
                            end
                        end
                    end
                end
                
                for user = 1:P.CDMAUsers
                    % We permute the himp to combine such that the H-1 will
                    % cancel out the interference inside rxsymbols
                    H = reshape(himp_perm(user,:,:,:),...
                        P.RakeFingers*P.RXperUser, P.TXperUser);
                    
                    % MIMO detector
                    switch P.Detector
                        case 'ZF'
                            % Zero Forcing Detector
                            desp_bits = pinv(H) * reshape(rxsymbols(user,:,:,:), P.RakeFingers*P.RXperUser, []);
                            
                            if strcmp(P.Mode,'HighDiversity')
                                desp_bits = sum(desp_bits,1); %summing for diversity
                            end
                            
                            switch P.Decision
                                case 'Hard'
                                    hard_desc_bits(user,:,:) = desp_bits < 0;
                                otherwise
                                    hard_desc_bits(user,:,:) = real(desp_bits);
                            end
                            
                        case 'MMSE'
                            % MMSE
                            Ps = 1 / (SNRlin);
                            G = (H' * H + (P.TXperUser / Ps) * eye(P.TXperUser)) \ H'; %N_T lower rate streams = TXperuser
                            desp_bits = G * reshape(rxsymbols(user,:,:,:), P.RakeFingers*P.RXperUser, []);
                            
                            if strcmp(P.Mode,'HighDiversity')
                                desp_bits = sum(desp_bits,1); %summing for diversity
                            end
                            
                            switch P.Decision
                                case 'Hard'
                                    hard_desc_bits(user,:,:)= desp_bits < 0;
                                otherwise
                                    hard_desc_bits(user,:,:)= real(desp_bits);
                            end
                            
                        case 'SIC'
                            H_k = H;
                            y_k = reshape(rxsymbols(user,:,:,:), P.RakeFingers*P.RXperUser, []);
                            for k = 1:P.TXperUser
                                
                                e = [1, zeros(1, P.TXperUser - k)];
                                g = e * pinv(H_k);
                                
                                % Estimate signal
                                s_hat = (g*y_k) < 0;
                                y_k = y_k - H(:, k)*s_hat;
                                
                                % Remove k-th column of H_k
                                H_k = H_k(:, 2:end);
                                
                                switch P.Decision %not good, to complete
                                    case 'Hard'
                                        hard_desc_bits(user,k,:) = s_hat;
                                    otherwise
                                end
                            end
                        otherwise
                            disp('Detector not supported')
                    end
                end
            otherwise
                disp('Source Encoding not supported')
        end

        if strcmp(P.Mode,'HighRate')
            decoded_bits = zeros(P.CDMAUsers,TX,NbTXBits/2);
            for user=1:P.CDMAUsers
                % UN-PN
                unPN_symbols = zeros(TX, NbTXBits);
                for t=1:TX
                    switch P.Decision
                        case 'Hard'
                            unPN_symbols(t,:) = xor(PNSequence(user,:), squeeze(hard_desc_bits(user,t,:)).');
                        otherwise
                            unPN_symbols(t,:) = (2*PNSequence(user,:)-1).* squeeze(hard_desc_bits(user,t,:)).';
                    end
                end
                
                % De-interleaver
                unPN_symbols=double(unPN_symbols);
                if strcmp(P.Interleaving, 'On')
                    unPN_symbols=matdeintrlv(unPN_symbols.',32,12).';
                end
                
                % Decoding Viterbi
                for t = 1:TX
                    decoded_bits(user,t,:) = convDec(unPN_symbols(t,:).').';
                end
            end
        else % Diversity mode
            for user=1:P.CDMAUsers
                % UN-PN
                switch P.Decision
                    case 'Hard'
                        unPN_symbols = xor(PNSequence(user,:), squeeze(hard_desc_bits(user,1,:)).');
                    otherwise
                         unPN_symbols = (2*PNSequence(user,:)-1) .* squeeze(hard_desc_bits(user,1,:)).';
                end
                % De-interleaver
                unPN_symbols=double(unPN_symbols);
                if strcmp(P.Interleaving, 'On')
                    unPN_symbols=matdeintrlv(unPN_symbols.',32,12).';
                end
                
                % Decoding Viterbi
                decoded_bits(user,1,:) = convDec(unPN_symbols.').';
            end
        end
        
        
        % Remove the 8 bit encoding trail
        for user=1:P.CDMAUsers
            rxbits(user,:,:) = decoded_bits(user,:,1:end-P.Q_Ind-8);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch P.Mode
            case 'HighRate'
                rxbits = reshape(rxbits, P.CDMAUsers, P.TXperUser, P.NumberOfBits);
            case 'HighDiversity'
                rxbits = reshape(rxbits, P.CDMAUsers, 1, P.NumberOfBits);
        end
        % Summing errors
        for user=1:P.CDMAUsers
            errors =  sum(sum(rxbits(user,:,:) ~= bits(user,:,:)));
            Results(ss) = Results(ss) + errors;        
        end
    end
end
if strcmp(P.Mode,'HighRate')
    BER = Results/(P.NumberOfBits*P.NumberOfFrames*TX*P.CDMAUsers);
else
    BER = Results/(P.NumberOfBits*P.NumberOfFrames*P.CDMAUsers);
end
end
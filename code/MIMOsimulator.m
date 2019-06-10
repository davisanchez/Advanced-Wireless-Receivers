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
    %NUsers = P.CDMAUsers;
    
    % Generate the spreading sequences
    HadamardMatrix = hadamard(P.HadLen)/sqrt(P.HadLen);    % Normalization      
    
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42); % 42 but could pick anything % TODO Magick numbers
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
                           
    PNSequence = step(LongCode).';
    
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

figure;
for frame = 1:P.NumberOfFrames
    
    frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bits = randi([0 1],TX,P.NumberOfBits); % Random Data, High rate mode
    
    % Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],TX,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = zeros(TX,NbTXBits);
    for t=1:TX
        encoded_bits(t,:) = convEnc(bits_ind(t,:).').';
    end
    
    % Symbol repetition
    % Symbols shall not be repeated for a data rate of 9600bps, the one
    % we are simulating (7-6)
    % encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    
    % Interleaver
    if strcmp(P.Interleaving, 'On')
        encoded_bits=matintrlv(encoded_bits.',32,12).';  % TODO Magick numbers
    end
    
    % Pulse Shape (PNSequence)
    PN_symbols = zeros(TX, NbTXBits);
    for t=1:TX
        PN_symbols(t,:) = xor(PNSequence, encoded_bits(t,:)); % TODO add user loop here
    end

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = zeros(TX,P.HadLen, NbTXBits);
    for t = 1:TX
        symbol_spread(t,:,:) = HadamSequence * symbols(t,:);
    end   
    waveform = reshape(symbol_spread,TX,NumberOfChips); 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = 1;
        case 'Multipath'
            himp = sqrt(1/2)* (randn(RX,TX,P.ChannelLength) +...
                            1i * randn(RX,TX,P.ChannelLength));
                        
            for r = 1:RX % Normalization for each combination
                for t = 1:TX
                    himp(r,t,:) = himp(r,t,:)/sqrt(sum(abs(himp(r,t,:)).^2)); 
                end
            end          
        otherwise
            disp('Channel not supported')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    
    switch P.ChannelType
        case {'AWGN','ByPass'}
            snoise = (randn(1,NumberOfChips) + ...
                      1i* randn(1,NumberOfChips) );
        case {'Multipath','Fading'}
            snoise = (randn(RX,NumberOfChipsRX) + ...
                      1i* randn(RX,NumberOfChipsRX) );    
    end
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass' 
                y = waveform * himp;
            case 'AWGN'
                y = waveform * himp + noise;  
            case 'Multipath'
                y = zeros(RX,NumberOfChipsRX); %Normally add the users here!
                for r = 1:RX
                    for t = 1:TX
                        y(r,:) = y(r,:) + conv(squeeze(waveform(t,:)),squeeze(himp(r,t,:)));  
                    end
                end
                y = y + noise;
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
                rxsymbols = zeros(RX,P.RakeFingers,NbTXBits); % TODO diversity or High rate possibility???
                
                % Separation between antennas ? How to get himp ?
                for r = 1:RX
                    for finger = 1:P.RakeFingers   
                        if strcmp(P.ChannelType,'Multipath')
                            % Despreading
                            rxsymbols(r,finger,:) = HadamSequence.'*reshape(...
                                y(r,finger:finger+NumberOfChips-1),...
                                P.HadLen, NbTXBits);
                        end
                    end
                end
            otherwise
                disp('Source Encoding not supported')
        end
        % We permute the himp to combine such that the H-1 will
        % cancel out the interference inside rxsymbols
        H = reshape(permute(himp(:,:,1:P.RakeFingers), [1,3,2]), P.RakeFingers*P.RXperUser, P.TXperUser);
                
        % MIMO detector
        switch P.Detector
            case 'ZF'
                % Zero Forcing Detector
                desp_bits = pinv(H) * reshape(rxsymbols, P.RakeFingers*P.RXperUser, []);
                % Hard Decision
                desc_bits = desp_bits < 0;
                
            case 'MMSE'
                % MMSE
                Ps = 1 / (SNRlin);
                G = (H' * H + (P.TXperUser / Ps) * eye(P.TXperUser)) \ H'; %N_T lower rate streams = TXperuser
                desp_bits = G * reshape(rxsymbols, P.RakeFingers*P.RXperUser, []);                
                % Hard Decision
                desc_bits = desp_bits < 0;
                
            case 'SIC'
                H_k = H;
                y_k = reshape(rxsymbols, P.RakeFingers*P.RXperUser, []);
                for k = 1:P.TXperUser
                    
                    e = [1, zeros(1, P.TXperUser - k)];
                    g = e * pinv(H_k);
                    
                    % Estimate signal
                    s_hat = (g*y_k) < 0;
                    y_k = y_k - H(:, k)*s_hat;
                    
                    % Remove k-th column of H_k
                    H_k = H_k(:, 2:end);
                    
                    %
                    desc_bits(k,:) = s_hat;
                end
            otherwise
                disp('Detector not supported')
        end

        % UN-PN
        unPN_symbols = zeros(TX, NbTXBits);
        for t=1:TX
            unPN_symbols(t,:) = xor(PNSequence, desc_bits(t,:)); % TODO add user loop here
        end
        
        % De-interleaver
        unPN_symbols=double(unPN_symbols);   
        if strcmp(P.Interleaving, 'On')
           unPN_symbols=matdeintrlv(unPN_symbols.',32,12).';
        end    
        
        % Decoding Viterbi
        decoded_bits = zeros(TX,NbTXBits/2);
        for t = 1:TX
            decoded_bits(t,:) = convDec(unPN_symbols(t,:).').'; % TODO, beurk beurk no??
        end
        
        % Remove the 8 bit encoding trail
        rxbits = decoded_bits(:,1:end-P.Q_Ind-8); % TODO magick number

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BER count
        errors =  sum(sum(rxbits ~= bits));
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(P.NumberOfBits*P.NumberOfFrames*TX);
end
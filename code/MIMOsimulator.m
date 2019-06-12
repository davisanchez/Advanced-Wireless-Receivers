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
    HadamardMatrix = hadamard(P.HadLen);     
    
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42); % 42 but could pick anything % TODO Magick numbers
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
     
    for user=1:P.CDMAUsers                       
        PNSequence(user,:) = step(LongCode).';
    end
    %change with each user
    
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

for frame = 1:P.NumberOfFrames
    
    frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch P.Mode
        case 'HighRate'
            bits = randi([0 1],TX,P.NumberOfBits); % Random Data, High rate mode  
            
            % Add Frame Quality Indicator (bonus)
            bits_ind = [bits randi([0 1],TX,P.Q_Ind)];
        case 'HighDiversity'
            bits = randi([0 1],1,P.NumberOfBits); % Few data, redundancy
            % Duplicate for each antenna
            bits_dupl = ones(TX,1)*bits;
            % Add Frame Quality Indicator (bonus)
            bits_ind = [bits_dupl randi([0 1],TX,P.Q_Ind)];
        otherwise
            disp('Mode not supported');
    end
    


    % Convolutional encoding
    encoded_bits = zeros(TX,NbTXBits, P.CDMAUsers);
    for t=1:TX
        for user=1:P.CDMAUsers
            encoded_bits(t,:,user) = convEnc(bits_ind(t,:).').';
        end
    end
    
    % Symbol repetition
    % Symbols shall not be repeated for a data rate of 9600bps, the one
    % we are simulating (7-6)
    % encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    
    % Interleaver
    if strcmp(P.Interleaving, 'On')
        for user=1:P.CDMAUsers
            encoded_bits(:,:,user)=matintrlv(encoded_bits(:,:,user).',P.IntrlvRows,NbTXBits/P.IntrlvRows).';
        end
    end
    
    % Pulse Shape (PNSequence)
    PN_symbols = zeros(TX, NbTXBits, P.CDMAUsers);
    for t=1:TX
        for user=1:P.CDMAUsers
            PN_symbols(t,:,user) = xor(PNSequence(user), encoded_bits(t,:,user)); % TODO add user loop here
        end
    end 

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1); 
    
    % Spreading with Hadamard
    symbol_spread = zeros(TX,P.HadLen, NbTXBits, P.CDMAUsers);
    for t = 1:TX
        for user=1:P.CDMAUsers
            symbol_spread(t,:,:,user) = HadamSequence * symbols(t,:,user);
        end
    end  
    %for user=1:P.CDMAUsers
        waveform = reshape(symbol_spread,TX,NumberOfChips,P.CDMAUsers); 
    %end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = 1;
        case 'Multipath'
            himp = sqrt(1/2)* (randn(RX,TX,P.ChannelLength) +...
                            1i * randn(RX,TX,P.ChannelLength));
        otherwise
            disp('Channel not supported')
    end
    %Normalization with Euclidean norm along third axis
    himp = himp./vecnorm(himp, 2, 3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    
    switch P.ChannelType
        case {'ByPass', 'AWGN'}
            snoise = (randn(1,NumberOfChips) + ...
                      1i* randn(1,NumberOfChips) );
        case 'Multipath'
            snoise = (randn(RX,TX,NumberOfChipsRX) + ...
                      1i* randn(RX,TX,NumberOfChipsRX) );    
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
                y = zeros(RX,NumberOfChipsRX, P.CDMAUsers); %Normally add the users here!
                for r = 1:RX
                    for t = 1:TX
                        for user=1:P.CDMAUsers
                            y(r,:,user) = y(r,:,user) + conv(squeeze(waveform(t,:,user)),squeeze(himp(r,t,:))) + squeeze(noise(r,t,:)).';
                        end
                    end
                end
                %y = y + noise;
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
                himp_perm = zeros(RX,P.RakeFingers,TX);
                               
                % Separation between antennas ? How to get himp ?
                for r = 1:RX
                    % Order the best combination of fingers (all contribution
                    % are combiend and unseparable !
                    [~,ind] = maxk(sum(abs(himp(r,:,:)),2),P.RakeFingers,3);
                    ind = squeeze(ind);
                    % Permute corresponding elements of himp for later
                    % calculation
                    himp_perm(r,:,:) = permute(himp(r,:,ind), [1,3,2]);
                    
                    for finger = 1:P.RakeFingers   
                        if strcmp(P.ChannelType,'Multipath')
                            % Despreading
                            rxsymbols(r,finger,:) = HadamSequence.'*reshape(...
                                y(r,ind(finger):ind(finger)+NumberOfChips-1),...
                                P.HadLen, NbTXBits);
                        end
                    end
                end
                
                % We permute the himp to combine such that the H-1 will
                % cancel out the interference inside rxsymbols
                H = reshape(himp_perm,...
                    P.RakeFingers*P.RXperUser, P.TXperUser);
                
                % MIMO detector
                switch P.Detector
                    case 'ZF'
                        % Zero Forcing Detector
                        desp_bits = pinv(H) * reshape(rxsymbols, P.RakeFingers*P.RXperUser, []);
                        
                        if strcmp(P.Mode,'HighDiversity')
                            desp_bits = sum(desp_bits,1); %summing for diversity
                        end
                        
                        % Hard Decision
                        desp_bits = desp_bits < 0;

                    case 'MMSE'
                        % MMSE
                        Ps = 1 / (SNRlin);
                        G = (H' * H + (P.TXperUser / Ps) * eye(P.TXperUser)) \ H'; %N_T lower rate streams = TXperuser
                        desp_bits = G * reshape(rxsymbols, P.RakeFingers*P.RXperUser, []);
                        
                        if strcmp(P.Mode,'HighDiversity')
                            desp_bits = sum(desp_bits,1); %summing for diversity
                        end
                        
                        % Hard Decision
                        desp_bits = desp_bits < 0;

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
                            desp_bits(k,:) = s_hat;
                        end
                    otherwise
                        disp('Detector not supported')
                end
            otherwise
                disp('Source Encoding not supported')
        end

        if strcmp(P.Mode,'HighRate')
            % UN-PN
            unPN_symbols = zeros(TX, NbTXBits, P.CDMAUsers);
            for t=1:TX
                for user=1:P.CDMAUsers
                    unPN_symbols(t,:, user) = xor(PNSequence(user,:), desp_bits(t,:)); % TODO add user loop here
                end
            end

            % De-interleaver
            unPN_symbols=double(unPN_symbols);   
            if strcmp(P.Interleaving, 'On')
                for user=1:P.CDMAUsers
                    unPN_symbols(:,:,user)=matdeintrlv(unPN_symbols(:,:,user).',32,12).';
                end
            end    

            % Decoding Viterbi
            decoded_bits = zeros(TX,NbTXBits/2,P.CDMAUsers);
            for t = 1:TX
                for user=1:P.CDMAUsers
                    decoded_bits(t,:,user) = convDec(unPN_symbols(t,:,user).').'; % TODO, beurk beurk no??
                end
            end
        else % Diversity mode
            % UN-PN
            for user=1:P.CDMAUsers
                unPN_symbols(:,:,user) = xor(PNSequence(:,user), desp_bits); % TODO add user loop here
                % De-interleaver
                unPN_symbols=double(unPN_symbols);   
                if strcmp(P.Interleaving, 'On')
                    unPN_symbols(:,:,user)=matdeintrlv(unPN_symbols(:,:,user).',32,12).';
                end    

                % Decoding Viterbi
                decoded_bits = convDec(unPN_symbols(:,:,user).').'; % TODO, beurk beurk no??
            end
        end
        
        
        % Remove the 8 bit encoding trail
        for user=1:P.CDMAUsers
            rxbits = decoded_bits(:,1:end-P.Q_Ind-8,user); % TODO magick number
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BER count
        errors =  sum(sum(rxbits ~= bits));
        Results(ss) = Results(ss) + errors;
        end
        
    end
end

BER = Results/(P.NumberOfBits*P.NumberOfFrames*TX*P.CDMAUsers); %NOT SURE
end
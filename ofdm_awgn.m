%{
Simulation of an OFDM wireless communication system for 4,16,64 QAM 
in the presence of AWGN channel
%}

clc;

%Number of subcarriers
N=4096;

%Length of cyclic prefix
cp_len=288;

%Various orders of modulation
mod=[4,16,64];

%Subcarrier spacing in Hz
scs=30e3;

%Total number of resources block
RB=273;

%Number of subcarriers per resource block
sc_rb=12;

%Total number of subcarriers
tot_sc=RB*sc_rb;

%Ratio of FFT size to the total number of utilised subcarriers 
ratio=ceil(N/tot_sc);

%Eb/No values
Eb_No=logspace(-1,3,100);

%Total number of iterations
Monte=5;

%Number of subcarriers per slot as per standards
sc_slot=14;

%List to store the BER values
ber=[];

for index=1:length(mod)
    
    %Modulation order
    mod_order=mod(index);

    %Number of bits per symbol
    k=log2(mod_order);

    %Total number of bits
    n_bits=k*tot_sc;

    %List to store SNR values 
    snr=[];

    for values=1:length(Eb_No)

        EbNo=Eb_No(values);

        snr_val = EbNo+10*log10(k)-10*log10(ratio);

        snr=[snr snr_val];

        ber_tot=0;

        ber_ratio=[];

        ber_num=[];

        for iterate=1:Monte

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transmitter side following operations takes place
%1. TX bits
%2. Modulation
%3. S2P
%4. IFFT
%5. P2S
%6. CP addition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Generating random bits
            bits=randi([0,1],n_bits,sc_slot);

            %List to store modulated symbols
            mod_data=[];

            for i=1:sc_slot

                %Modulating the input symbols
                modulated_1=qammod(bits(:,i),mod_order,"gray","InputType",...
                "bit","UnitAveragePower",false);

                mod_data=[mod_data modulated_1];

            end

            %Filling the entire bandwidth with zeros
            sc_pos=zeros(sc_slot,N);

            %Taking transpose of the modulated matrice
            modulated_symbols_tr=mod_data';

            %Only a fraction of the entire bandwidth will have data
            %carrying subcarriers rest will be guard band
            sc_pos(:,1:tot_sc)=modulated_symbols_tr;

            %S2P
            sc_prl=sc_pos';

            %Taking IFFT after serial to parallel conversion
            ifft_out=ifft(sc_prl,N);

            %P2S
           srl_vec_tx=ifft_out(:);

           srl_vec_tx=srl_vec_tx';

           %CP addition
           srl_vec_tx_cp=[];

           for i=1:sc_slot

               cp_add=srl_vec_tx(N*i-(cp_len-1):N*i);
               srl_vec_tx_cp=[srl_vec_tx_cp cp_add srl_vec_tx(N*(i-1)+1:N*i)];

           end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Channel- ideal wireless channel is being assumed and AWGN noise is added
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %adding AWGN noise to the transmitted symbols
           chn_awgn=awgn(srl_vec_tx_cp,snr_val,"measured");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Receiver side following operations takes place
%1. CP Removal
%2. S2P
%3. FFT
%4. P2S
%5. Demodulation
%6. RX bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           %CP removal
           row_vec_rx_wocp=[];

           for i=1:sc_slot
               var_1=(N+cp_len)+(i-1)*(N+cp_len);
               row_vec_rx_wocp=[row_vec_rx_wocp chn_awgn(var_1-N+1:var_1)];
           end

           %S2P
           row_vec_rx_prl=row_vec_rx_wocp';

           %Converting parallel data to matrice form
           mat_rx_data=reshape(row_vec_rx_prl,N,sc_slot);

           %Taking FFT of the resultant matrice
           fft_out=fft(mat_rx_data,N);

           %P2S
           fft_out_srl=fft_out';

           %Mapping to resource elements
           rx_re_data=fft_out_srl(:,1:tot_sc);

           %S2P
           rx_re_data_prl=rx_re_data';

           %Demodulation

           %List to store demodulated received bits
           received_bits=[];

           for i=1:sc_slot
               demodulated_1=qamdemod(rx_re_data_prl(:,i),mod_order,"gray",...
                   "OutputType","bit");
               received_bits=[received_bits demodulated_1];
           end

           %Converting the received bits to column form
           received_bits=received_bits(:);

           %Converting the input bits to column form
           input_bits=bits(:);

           %Calculating the BER
           [ber_1 ber_2]=biterr(input_bits,received_bits);

           ber_tot=ber_tot+ber_2;
            
        end

        %Averaging BER
        ber=[ber ber_tot/Monte];

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the theoretical BER plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theoretical_ber=[];

for i=1:length(mod)

    modulation=mod(i);
    y = qfunc(sqrt(3*log2(modulation)/(modulation-1)*Eb_No));
    y = 4*(1 - 1/sqrt(modulation))*y;
    y = y/log2(modulation);
    theoretical_ber = [theoretical_ber, y];

end
theoretical_ber_1=theoretical_ber(1:length(Eb_No));

theoretical_ber_2=theoretical_ber(length(Eb_No)+1:2*length(Eb_No));

theoretical_ber_3=theoretical_ber(2*length(Eb_No)+1:3*length(Eb_No));

semilogy(10*log10(Eb_No),theoretical_ber_1,'b--o');
hold on;
semilogy(10*log10(Eb_No),theoretical_ber_2,'g--o');
hold on;
semilogy(10*log10(Eb_No),theoretical_ber_3,'m--o');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating the simulation BER plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ber_1=ber(1:length(Eb_No));

ber_2=ber(length(Eb_No)+1:2*length(Eb_No));

ber_3=ber(2*length(Eb_No)+1:3*length(Eb_No));

%Plot
semilogy(10*log10(Eb_No),ber_1,'b--^');
hold on;
semilogy(10*log10(Eb_No),ber_2,'g--^');
hold on;
semilogy(10*log10(Eb_No),ber_3,'m--^');
ylim([1e-6, 1e1]);
xlabel("$10\log\left(\frac{E_b}{N_0}\right)$","Interpreter", "latex");
ylabel("BER");
title("BER waterfall curves")
legend("4 QAM theory", "16 QAM theory", "64 QAM theory",...
    "4 QAM sim", "16 QAM sim", "64 QAM sim",...
    "NumColumns", 2);
grid on;






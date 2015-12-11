function BBQ=slideFFT(BBQ,start_turns,end_turns,shift_turns,window_turns,method,nharm)
 if ~exist('nharm','var')
     % third parameter does not exist, so default it to something
      nharm = 1;
 end
 
BBQ_fft=BBQ;
BBQ_fft.time=BBQ_fft.time(start_turns:end_turns);
BBQ_fft.turns=BBQ_fft.turns(start_turns:end_turns);
BBQ_fft.data=BBQ_fft.data(start_turns:end_turns);
L = window_turns; % width of window
N = length(BBQ_fft.data); % number of windows
del=mod(N,L);
BBQ_fft.turns=BBQ_fft.turns(1:end-del);
BBQ_fft.data=BBQ_fft.data(1:end-del);
N = length(BBQ_fft.data); 
BBQ_fft.slide_turns=1:shift_turns:N-L;
BBQ_fft.slide_time=BBQ_fft.time(BBQ_fft.slide_turns);

freq_vec=[];
amp_vec=[];
tune_vec=[];

for ii=1:shift_turns:N-L
    BBQ_fft.data_table=BBQ_fft.data(ii:ii+L-1)';

    if strcmp(method,'NAFF')
        disp(['performing NAFF on ',num2str(ii),'/',num2str(N-L),' windows of ',num2str(L),' turns']);
        [~,~,~,~,~,~,frespec,~]=freqan(BBQ_fft.data_table,nharm,'all');

        freq=reshape(frespec(1:end-nharm,3)',nharm, length(frespec(1:end-nharm,3))/nharm);
        freq=freq.*sign(freq);
        amp=reshape(frespec(1:end-nharm,4)',nharm, length(frespec(1:end-nharm,4))/nharm);

        amp(freq<0.2)=nan; % clean low frequency stuff
        freq(freq<0.2)=nan;

        [~,ind]=max(amp,[],1);
        tune=freq(ind);
        


    elseif strcmp(method,'FFT')
        disp(['performing FFT on ',num2str(ii),'/',num2str(N-L),' windows of ',num2str(L),' turns']);
        specFFT=fftshift(fft(BBQ_fft.data_table));
        amp=abs(specFFT);
        freq=linspace(-1/2,1/2,L);
        tune=[];
    %     ind=find(freq>0.2);
    %     amp_pos=amp(ind);
    %     freq_pos=freq(ind);
    %     [~,indmax]=max(amp_pos);
    %     tune=freq_pos((indmax));
    %     disp(['tune ',num2str(tune)])
    %     if isempty(indmax); tune=nan; end
    %     tune_vec=[tune_vec,tune];
    end

    freq_vec=[freq_vec,freq(:)];
    amp_vec=[amp_vec,amp(:)];
    tune_vec=[tune_vec,tune];
end

if strcmp(method,'FFT')
    BBQ.FFT_slide_turns=BBQ_fft.slide_turns;
    BBQ.FFT_slide_time=BBQ_fft.slide_time;
    BBQ.FFT_slide_data=BBQ_fft.data(BBQ_fft.slide_turns);
    BBQ.FFT_slide_tune=tune_vec;
    BBQ.FFT_freq=freq_vec;
    BBQ.FFT_amp=amp_vec;
else
    BBQ.NAFF_slide_turns=BBQ_fft.slide_turns;
    BBQ.NAFF_slide_time=BBQ_fft.slide_time;
    BBQ.NAFF_slide_data=BBQ_fft.data(BBQ_fft.slide_turns);
    BBQ.NAFF_slide_tune=tune_vec;
    BBQ.NAFF_freq=freq_vec;
    BBQ.NAFF_amp=amp_vec;
end
end



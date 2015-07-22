function BBQ=slideFFT(BBQ,start_turns,end_turns,window_turns,method,nharm)
 if ~exist('nharm','var')
     % third parameter does not exist, so default it to something
      nharm = 1;
 end
 
BBQ_fft=BBQ;

BBQ_fft.turns=BBQ_fft.turns(start_turns:end_turns);
BBQ_fft.data=BBQ_fft.data(start_turns:end_turns);

L = window_turns; % width of window
N = length(BBQ_fft.data); % number of windows
mi=repmat(1:L,[N-L+1,1])+repmat((0:(N-L))',[1,L]);
BBQ_fft.data_table = BBQ_fft.data(mi);
BBQ_fft.slide_turns=1:size(BBQ_fft.data_table,1);


if strcmp(method,'NAFF')
    
    [~,~,~,~,~,~,frespec,~]=freqan(BBQ_fft.data_table',nharm,'all');
    
    freq=reshape(frespec(1:end-nharm,3)',nharm, length(frespec(1:end-nharm,3))/nharm );
    freq=freq.*sign(freq);
    amp=reshape(frespec(1:end-nharm,4)',nharm, length(frespec(1:end-nharm,4))/nharm );
   
    amp(freq<0.2)=nan; % clean low frequency stuff
    freq(freq<0.2)=nan;
    
    [~,ind]=max(amp,[],1);
    tune_vec=[];
    for ii=1:length(freq)
        tune_vec(ii)=freq(ind(ii),ii);
    end


elseif strcmp(method,'FFT')
    disp(['performing FFT on ',num2str(size(BBQ_fft.data_table,1)),' windows of ',num2str(size(BBQ_fft.data_table,2)),' turns']);
    specFFT=fftshift(fft(BBQ_fft.data_table,[],2));
    amp=abs(specFFT);
    freq=linspace(-1/2,1/2,L);
    tune_vec=[];
%     ind=find(freq>0.2);
%     amp_pos=amp(ind);
%     freq_pos=freq(ind);
%     [~,indmax]=max(amp_pos);
%     tune=freq_pos((indmax));
%     disp(['tune ',num2str(tune)])
%     if isempty(indmax); tune=nan; end
%     tune_vec=[tune_vec,tune];
end

BBQ.slide_turns=BBQ_fft.slide_turns;
BBQ.slide_data=BBQ_fft.data(BBQ_fft.slide_turns);
BBQ.slide_tune=tune_vec;
BBQ.freq=freq;
BBQ.amp=amp;

end



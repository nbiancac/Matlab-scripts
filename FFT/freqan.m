function [avetune,errtune,tunes,phaseadv,amplit,numturns,frespec,nkick] = freqan(data,nterms,cas)
% FREQAN    : Refined Fourier analysis from turn-by-turn BPM data
% Usage     :  [avetune,errtune,tunes,phaseadv,amplit,numturns,frespec,nkick] = freqan(data,nterms,cas)
% data      : BPM data matrix (each column a different BPM)
% nterm     : scalar indicating the number of terms
% cas       : string for indicating type of analysis 
%             i)  cas='all', then analyze all data
%             ii) cas anything, else then analyze up to decoherence point
% avetune   : average tune for all BPMs
% errtune   : statistical error in tune estimation
% tunes     : array with all the tunes for each BPM
% phaseadv  : array with the phase advance for each BPM
% amplit    : array with amplitude associated to the highest pick
%             for each BPM
% numturns  : number of analyzed turns

		
% test the number and the class of the arguments
switch nargin
 case 0
  disp ('Please, give me some data to analyse!');
  exit
 case 1
  nterms=1;
  cas=0;
 case 2
  if isnumeric(nterms);
    cas=0;
    if nterms==0
      cas=nterms;
      nterms=1;
    end
  else
    cas=nterms;
    nterms=1;
  end
 otherwise
   if nterms==0
     nterms=1;
   end
end

% check for bad BPM readings 
anynoisy=any(isnan(data));                 
allnoisy=all(isnan(data));
goodbpms=all(isfinite(data));

% calculate total number of turns and BPMs 
[allnumturns,allnumbpm]=size(data);        
shft=[1 1:allnumturns-1];                                                   
lastbpm=max(find(all(isfinite(data))));
dataend=data(shft,lastbpm);	 

% get rid of "noisy" BPMs 
data=[dataend data(:,goodbpms)];          
% data=[data(:,goodbpms)];          

% get strings with number of BPMs  
numbpm=size(data,2)-1;
 
% find the turn where the kick occurs 
if cas~='all'
  [maxrms,nkick]=max(abs(diff(std(data'))));
  if nkick < 1000 
    % get rid of the "noisy" data                           
    data(1:nkick,:)=[];	
    disp(nkick)
  else
    disp(['nkick = ' num2str(nkick) ' undefined.  Analyzing all data!']);
    nkick=1;
  end 
    % find the coherence
  threshold=0.1*maxrms;                      
  on=find(std(data')> threshold);
  coher=on(length(on));
  
  % get rid of the incoherent turns
  if coher  < allnumturns-nkick            
     data(coher+1:allnumturns-nkick,:)=[];
  end 
else
  nkick=1;
  coher=size(data,1);
end
  
% get strings with number of turns and BPMs 
numturns=size(data,1); 

disp(['Frequency analysis of ' num2str(size(data,1)) ' turns from ' ...
      num2str(numbpm) ' data sets' ]); 

% call the frequency analysis "matlab" interface function 

result=matnafterms(data,nterms);   


% !!!!!! results !!!!!!
mytunes=NaN*ones(1,allnumbpm);               
myphases=NaN*ones(1,allnumbpm);             
myamplit=NaN*ones(1,allnumbpm);


mytunes(goodbpms)=result(1+nterms:nterms:end,3);
myphases(goodbpms)=cumsum(mod(diff(sign(result(1:nterms:end,3)).* ...
				   result(1:nterms:end,5))/2/pi,1/2)); 
myamplit(goodbpms)=result(1+nterms:nterms:end,4);

frespec=NaN*ones(allnumbpm*nterms,7);
bpmsfrespec=reshape(1:allnumbpm*nterms,nterms,allnumbpm);
goodbpmsterms=reshape(bpmsfrespec(:,goodbpms),nterms*sum(goodbpms), ...
		      1);
frespec(goodbpmsterms,:)=result(1+nterms:end,:);
% frespec=result;

% check the tune error and proceed accordingly 
errtune=std(abs(mytunes(goodbpms)));
if errtune < 1/sqrt(numbpm)/numturns       
   avetune=mean(abs(mytunes(goodbpms)));
   goodtunes=goodbpms;
else
   disp(['Error in the tune calculation too large, recomputing' ...
	     ' average']); 
   disp(['Check for bad BPM readings or coupling']);
   goodtunes=(abs(abs(mytunes)-mean(abs(mytunes(goodbpms))))< errtune);
   errtune=std(abs(mytunes(goodtunes)));
   avetune=mean(abs(mytunes(goodtunes)));
end 
tunes=NaN*ones(1,allnumbpm);
phaseadv=NaN*ones(1,allnumbpm);
amplit=NaN*ones(1,allnumbpm);
tunes(goodbpms)=mytunes(goodbpms);
phaseadv(goodtunes)=myphases(goodtunes);
amplit(goodtunes)=myamplit(goodtunes);
disp(['Average tune = ' num2str(avetune) ', error =  ' num2str(errtune)]);
frespec=result;


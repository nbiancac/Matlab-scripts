function L=unwrap_ADT(L)

L.data=[];
L.bucket=[];
for ii=1:size(L.raw,1)
    if sum(L.raw(ii,:))==0
        
    else
        data=L.raw(ii,:);
        data(data==0)=nan;
        L.data=[L.data;double(data(1:end))];
        L.bucket=[L.bucket,ii];
    end
    
end



end
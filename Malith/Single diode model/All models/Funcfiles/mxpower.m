 function [Vm,Im,max_index] = mxpower(Voc_index,Isc_index,I,V)
    %Find Im and Vm from MPP tracking
    %Calculating power for the reverse bias and forward bias measurements
    %may become an issue because indexes are swapped. e.g in forward I(1)
    %will contain values near Isc. But in reverse bias I(1) will be
    %negative or closer to zero
    %First if statement is for forward bias and the else is for negative
    if (Voc_index > Isc_index )
        mapvect = Isc_index:1:Voc_index;
        mpp = abs(I(mapvect).*V(mapvect));
        max_index = find(mpp==max(mpp)) + Isc_index;
        Im = I(max_index);
        Vm = V(max_index);
    
    else 
        mapvect = Voc_index:1:Isc_index;
        mpp = abs(I(mapvect).*V(mapvect));
        max_index = find(mpp==max(mpp)) + Voc_index;
        Im = I(max_index);
        Vm = V(max_index);

    end
    
    %find gradient at MPP; should be close to zero
        


end 
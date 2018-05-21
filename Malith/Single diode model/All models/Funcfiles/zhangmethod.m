function [dvdi] = zhangmethod(V,I,mapvect,Voc_index,Isc_index)
   %assumption works for active quarant
    if (Voc_index > Isc_index )
    
    c = Voc_index;
    
    c1 = (V(c-3) - V(c-4))/(I(c-3) - I(c-4));
    c2 = (V(c-2) - V(c-3))/(I(c-2) - I(c-3));
    c3 = (V(c-1) - V(c-2))/(I(c-1) - I(c-2));
    c4 = (V(c)   - V(c-1))/(I(c) - I(c-1));
    c5 = (V(c+1) - V(c))/(I(c+1) - I(c));
    c6 = (V(c+2) - V(c+1))/(I(c+2) - I(c+1));
    
    dvdi = [c1,c2,c3,c4,c5,c6];
    
    
    
    %polyfit()
    
    else 
        
    c = mapvect;
    c1 = -(V(c-3) - V(c-4))/(I(c-3) - I(c-4));
    c2 = -(V(c-2) - V(c-3))/(I(c-2) - I(c-3));
    c3 = -(V(c-1) - V(c-2))/(I(c-1) - I(c-2));
    c4 = -(V(c)   - V(c-1))/(I(c) - I(c-1));
    c5 = -(V(c+1) - V(c))/(I(c+1) - I(c));
    c6 = -(V(c+2) - V(c+1))/(I(c+2) - I(c+1));
    dvdi = [c1,c2,c3,c4,c5,c6];

    end
    



end
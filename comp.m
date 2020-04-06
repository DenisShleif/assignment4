classdef comp
    properties
        name
        n1
        n2
        val
        num
        rel
        subCircNum
    end
    
    methods
        function obj = comp(Name,N1,N2,Val,Num,Rel)
            obj.name = Name;
            obj.n1 = N1;
            obj.n2 = N2;
            obj.val = Val;
            obj.num = Num;
            if nargin == 6
               obj.rel = Rel; 
            end
        end
        
        function [out] = getUnit(obj)
            out = obj.name;
            switch obj.getType()
                case 'v'
                    out = [out,' (V)'];
                case 'i'
                    out = [out,' (A)'];
                case 'r'
                    out = [out,' (R)'];
                case 'l'
                    out = [out,' (H)'];
                case 'c'
                    out = [out,' (F)'];
            end
        end
        
        function [t] = getType(obj)
           t = obj.name(1);
        end
    end
end

